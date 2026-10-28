#include "dqmc/dqmc.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace dqmc {
namespace {

const char kMagic[8] = {'D', 'Q', 'M', 'C', 'C', 'K', 'P', 'T'};

std::uint64_t fnv1a_update(std::uint64_t hash, const void* data,
                           std::size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

template <typename T>
void put(std::string& buffer, std::uint64_t& hash, const T& value) {
  buffer.append(reinterpret_cast<const char*>(&value), sizeof(T));
  hash = fnv1a_update(hash, &value, sizeof(T));
}

void put_blob(std::string& buffer, std::uint64_t& hash, const void* data,
              std::uint64_t size) {
  buffer.append(static_cast<const char*>(data), size);
  hash = fnv1a_update(hash, data, size);
}

class Reader {
 public:
  Reader(const char* data, std::size_t size) : data_(data), size_(size) {}

  template <typename T>
  T get() {
    T value;
    if (offset_ + sizeof(T) > size_)
      throw std::runtime_error("checkpoint truncated");
    std::memcpy(&value, data_ + offset_, sizeof(T));
    hash_ = fnv1a_update(hash_, &value, sizeof(T));
    offset_ += sizeof(T);
    return value;
  }

  void get_blob(void* dest, std::uint64_t size) {
    if (offset_ + size > size_) throw std::runtime_error("checkpoint truncated");
    std::memcpy(dest, data_ + offset_, size);
    hash_ = fnv1a_update(hash_, dest, size);
    offset_ += size;
  }

  std::uint64_t hash() const { return hash_; }
  std::size_t offset() const { return offset_; }

 private:
  const char* data_;
  std::size_t size_;
  std::size_t offset_ = 0;
  std::uint64_t hash_ = 14695981039346656037ULL;
};

}  // namespace

void Solver::save_checkpoint(const std::string& path,
                             std::uint64_t sweep_tag) const {
  if (!params_.recompute_each_slice && !params_.recompute_after_sweep)
    throw std::runtime_error(
        "checkpointing requires a recompute flag so the boundary Green "
        "function is reproducible");

  std::string buffer;
  std::uint64_t hash = 14695981039346656037ULL;
  put_blob(buffer, hash, kMagic, sizeof(kMagic));
  put(buffer, hash, checkpoint_version);
  put(buffer, hash, params_.lx);
  put(buffer, hash, params_.ly);
  put(buffer, hash, params_.slices);
  put(buffer, hash, params_.t);
  put(buffer, hash, params_.abs_u);
  put(buffer, hash, params_.mu);
  put(buffer, hash, params_.dtau);
  put(buffer, hash, params_.seed);
  std::uint8_t flags = (params_.recompute_each_slice ? 1 : 0) |
                       (params_.recompute_after_sweep ? 2 : 0);
  put(buffer, hash, flags);
  put(buffer, hash, accepted_);
  put(buffer, hash, attempted_);
  put(buffer, hash, sweep_tag);

  const std::uint64_t field_count = fields_.size();
  put(buffer, hash, field_count);
  put_blob(buffer, hash, fields_.data(), field_count);

  for (const Matrix& b : bmat_)
    put_blob(buffer, hash, b.data(), sizeof(double) * b.size());

  std::ostringstream rng_stream;
  rng_stream << rng_;
  const std::string rng_text = rng_stream.str();
  const std::uint64_t rng_size = rng_text.size();
  put(buffer, hash, rng_size);
  put_blob(buffer, hash, rng_text.data(), rng_size);

  buffer.append(reinterpret_cast<const char*>(&hash), sizeof(hash));

  const std::string tmp_path = path + ".tmp";
  {
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot open checkpoint for writing");
    out.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    if (!out) throw std::runtime_error("checkpoint write failed");
  }
  if (std::rename(tmp_path.c_str(), path.c_str()) != 0)
    throw std::runtime_error("checkpoint rename failed");
}

Solver Solver::load_checkpoint(const std::string& path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("cannot open checkpoint");
  const std::streamsize size = in.tellg();
  in.seekg(0);
  std::string buffer(static_cast<std::size_t>(size), '\0');
  in.read(buffer.data(), size);
  if (!in) throw std::runtime_error("checkpoint read failed");

  if (buffer.size() < sizeof(kMagic) + sizeof(std::uint64_t))
    throw std::runtime_error("checkpoint too small");
  Reader reader(buffer.data(), buffer.size());

  char magic[8];
  reader.get_blob(magic, sizeof(magic));
  if (std::memcmp(magic, kMagic, sizeof(kMagic)) != 0)
    throw std::runtime_error("not a DQMC checkpoint");
  const std::uint32_t version = reader.get<std::uint32_t>();
  if (version != checkpoint_version)
    throw std::runtime_error("unsupported checkpoint version");

  Parameters params;
  params.lx = reader.get<int>();
  params.ly = reader.get<int>();
  params.slices = reader.get<int>();
  params.t = reader.get<double>();
  params.abs_u = reader.get<double>();
  params.mu = reader.get<double>();
  params.dtau = reader.get<double>();
  params.seed = reader.get<std::uint64_t>();
  const std::uint8_t flags = reader.get<std::uint8_t>();
  params.recompute_each_slice = (flags & 1) != 0;
  params.recompute_after_sweep = (flags & 2) != 0;

  Solver solver(params, UninitializedTag{});
  solver.accepted_ = reader.get<std::uint64_t>();
  solver.attempted_ = reader.get<std::uint64_t>();
  solver.sweep_tag_ = reader.get<std::uint64_t>();

  const std::uint64_t field_count = reader.get<std::uint64_t>();
  if (field_count !=
      static_cast<std::uint64_t>(solver.nsites_) * params.slices)
    throw std::runtime_error("checkpoint field count mismatch");
  solver.fields_.resize(field_count);
  reader.get_blob(solver.fields_.data(), field_count);
  for (const auto field : solver.fields_)
    if (field != -1 && field != 1)
      throw std::runtime_error("checkpoint field value out of range");

  solver.bmat_.assign(params.slices, Matrix(solver.nsites_, solver.nsites_));
  for (Matrix& b : solver.bmat_)
    reader.get_blob(b.data(), sizeof(double) * b.size());

  const std::uint64_t rng_size = reader.get<std::uint64_t>();
  std::string rng_text(rng_size, '\0');
  reader.get_blob(rng_text.data(), rng_size);

  std::uint64_t stored_hash;
  if (reader.offset() + sizeof(stored_hash) != buffer.size())
    throw std::runtime_error("checkpoint size mismatch");
  std::memcpy(&stored_hash, buffer.data() + reader.offset(),
              sizeof(stored_hash));
  if (stored_hash != reader.hash())
    throw std::runtime_error("checkpoint checksum mismatch");

  std::istringstream rng_stream(rng_text);
  rng_stream >> solver.rng_;
  if (!rng_stream) throw std::runtime_error("checkpoint RNG state unreadable");

  solver.green_ = solver.green_at_slice(0);
  return solver;
}

}  // namespace dqmc
