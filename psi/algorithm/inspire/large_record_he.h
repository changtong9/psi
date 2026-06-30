#pragma once

#include <cstdint>
#include <vector>

#include "hexl/hexl.hpp"

namespace psi::inspire::lrhe {

constexpr uint64_t kLrCrtQ1 = 268369921ULL;
constexpr uint64_t kLrCrtQ2 = 249561089ULL;
constexpr uint64_t kLrCrtMod = kLrCrtQ1 * kLrCrtQ2;
constexpr uint64_t kLrRootOfUnityCrt4096 = 3375402822066082ULL;
constexpr uint64_t kLrRootOfUnityCrt2048 = 38878761190133527ULL;

class Poly {
 public:
  Poly();
  Poly(uint64_t len, uint64_t modulus, bool ntt = false);
  Poly(const std::vector<uint64_t>& vec, uint64_t modulus, bool ntt);

  uint64_t get_length() const;
  uint64_t get_modulus() const;
  bool get_nttform() const;
  void set_nttform(bool ntt);
  uint64_t* data();
  const uint64_t* data() const;
  uint64_t& operator[](size_t idx);
  const uint64_t& operator[](size_t idx) const;
  std::vector<uint64_t>& get_data();
  const std::vector<uint64_t>& get_data() const;

 private:
  uint64_t length_;
  uint64_t modulus_;
  bool nttform_;

 public:
  std::vector<uint64_t> payload;
};

void PolyMult(Poly& lhs, Poly& rhs, Poly& result);
void PolyAdd(Poly& lhs, Poly& rhs, Poly& result);
void PolySub(Poly& lhs, Poly& rhs, Poly& result);
void PolyFMA(Poly& mul, uint64_t constant, Poly& result);
void PolyToNTT(Poly& poly, intel::hexl::NTT& ntts);
void PolyToCoef(Poly& poly, intel::hexl::NTT& ntts);

class Secret {
 public:
  Secret(std::vector<int64_t>& vec, uint64_t modulus, bool ntt);
  intel::hexl::NTT& get_ntt();
  Poly& get_data();
  const Poly& get_data() const;
  uint64_t get_modulus() const;
  uint64_t get_length() const;
  bool get_nttform() const;
  void to_ntt_form();

 private:
  intel::hexl::NTT ntts_;
  Poly data_;
};

class RlweCiphertext {
 public:
  RlweCiphertext();
  RlweCiphertext(uint64_t len, uint64_t modulus, bool ntt);

  uint64_t get_degree() const;
  uint64_t get_modulus() const;
  void set_nttform(bool ntt);

  Poly a;
  Poly b;
};

class Plaintext {
 public:
  Plaintext(uint64_t len, uint64_t pmod);
  Plaintext(const std::vector<uint64_t>& vec, uint64_t pmod, bool ntt);

  uint64_t get_pmod() const;
  bool get_isntt() const;

  Poly pt;
};

class ApproximateRgswCiphertext {
 public:
  ApproximateRgswCiphertext(uint64_t len, uint64_t modulus, uint64_t b,
                            uint64_t z, uint64_t t);

  uint64_t get_degree() const;
  uint64_t get_modulus() const;
  uint64_t get_b() const;
  uint64_t get_z() const;
  uint64_t get_t() const;

 private:
  uint64_t degree_;
  uint64_t modulus_;
  uint64_t b_;
  uint64_t z_;
  uint64_t t_;

 public:
  std::vector<RlweCiphertext> ct_m;
  std::vector<RlweCiphertext> ct_sm;
};

void RlweBfvEncrypt(Plaintext& pt, RlweCiphertext& ct, Secret& sk, float sig);
void RlweBfvEncode(Plaintext& pt, RlweCiphertext& ct, Secret& sk, float sig);
void RlweBfvDecrypt(RlweCiphertext& ct, Plaintext& pt, Secret& sk);
void RlweBfvAdd(RlweCiphertext& add1, RlweCiphertext& add2,
                RlweCiphertext& result);
void ExternalProduct(const ApproximateRgswCiphertext& ct_rgsw,
                     RlweCiphertext& ct_bfv, RlweCiphertext& result,
                     intel::hexl::NTT& ntts);
void RgswEncode(Plaintext& pt, ApproximateRgswCiphertext& ct, Secret sk,
                float sig);

}  // namespace psi::inspire::lrhe
