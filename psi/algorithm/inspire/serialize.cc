#include "psi/algorithm/inspire/serialize.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "yacl/base/exception.h"

namespace psi::inspire {
namespace {

template <typename T>
void AppendPod(std::string& out, const T& value) {
  out.append(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
T ReadPod(const char*& ptr, const char* end) {
  YACL_ENFORCE_GE(end - ptr, static_cast<ptrdiff_t>(sizeof(T)));
  T value;
  std::memcpy(&value, ptr, sizeof(T));
  ptr += sizeof(T);
  return value;
}

template <typename T>
void AppendVector(std::string& out, const std::vector<T>& values) {
  AppendPod<uint64_t>(out, values.size());
  if (!values.empty()) {
    out.append(reinterpret_cast<const char*>(values.data()),
               values.size() * sizeof(T));
  }
}

template <typename T>
std::vector<T> ReadVector(const char*& ptr, const char* end) {
  const auto size = ReadPod<uint64_t>(ptr, end);
  YACL_ENFORCE_GE(end - ptr, static_cast<ptrdiff_t>(size * sizeof(T)));
  std::vector<T> out(size);
  if (size > 0) {
    std::memcpy(out.data(), ptr, size * sizeof(T));
    ptr += size * sizeof(T);
  }
  return out;
}

template <typename T>
void AppendVector2D(std::string& out,
                    const std::vector<std::vector<T>>& values) {
  AppendPod<uint64_t>(out, values.size());
  for (const auto& inner : values) {
    AppendVector(out, inner);
  }
}

template <typename T>
std::vector<std::vector<T>> ReadVector2D(const char*& ptr, const char* end) {
  const auto outer = ReadPod<uint64_t>(ptr, end);
  std::vector<std::vector<T>> out;
  out.reserve(outer);
  for (uint64_t i = 0; i < outer; ++i) {
    out.push_back(ReadVector<T>(ptr, end));
  }
  return out;
}

yacl::Buffer ToBuffer(const std::string& bytes) {
  return yacl::Buffer(bytes.data(), bytes.size());
}

}  // namespace

yacl::Buffer SerializeQuery(const InspireQuery& query) {
  std::string bytes;
  AppendVector(bytes, query.qu0);
  AppendVector(bytes, query.qu1);
  AppendVector2D(bytes, query.ksk_b);
  return ToBuffer(bytes);
}

InspireQuery DeserializeQuery(const yacl::ByteContainerView& buffer) {
  const char* ptr = reinterpret_cast<const char*>(buffer.data());
  const char* end = ptr + buffer.size();

  InspireQuery query;
  query.qu0 = ReadVector<uint64_t>(ptr, end);
  query.qu1 = ReadVector<uint64_t>(ptr, end);
  query.ksk_b = ReadVector2D<uint64_t>(ptr, end);
  YACL_ENFORCE(ptr == end, "unexpected trailing bytes in InspireQuery");
  return query;
}

yacl::Buffer SerializeResponse(const InspireResponse& response) {
  std::string bytes;
  AppendVector2D(bytes, response.doublepir_response);
  return ToBuffer(bytes);
}

InspireResponse DeserializeResponse(const yacl::ByteContainerView& buffer) {
  const char* ptr = reinterpret_cast<const char*>(buffer.data());
  const char* end = ptr + buffer.size();

  InspireResponse response;
  response.doublepir_response = ReadVector2D<uint64_t>(ptr, end);
  YACL_ENFORCE(ptr == end, "unexpected trailing bytes in InspireResponse");
  return response;
}

}  // namespace psi::inspire
