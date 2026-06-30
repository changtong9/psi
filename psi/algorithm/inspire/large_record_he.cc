#include "psi/algorithm/inspire/large_record_he.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

namespace psi::inspire::lrhe {

namespace {

void SampleRandom(Poly& a, uint64_t modulus, uint64_t length) {
  for (uint64_t i = 0; i < length; ++i) {
    uint64_t value = (static_cast<uint64_t>(std::rand()) << 31) |
                     static_cast<uint64_t>(std::rand());
    a[i] = value % modulus;
  }
}

void SampleGauss(std::vector<uint64_t>& err, float st_dev, uint64_t modulus) {
  std::default_random_engine engine(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<double> gaussian(0.0, st_dev);
  for (size_t i = 0; i < err.size(); ++i) {
    int64_t value = static_cast<int64_t>(std::llround(gaussian(engine)));
    err[i] = (value < 0) ? (modulus - static_cast<uint64_t>(-value))
                         : static_cast<uint64_t>(value);
  }
}

void ApproximateGadgetDecomp(Poly& poly, std::vector<Poly>& mat, uint64_t b,
                             uint64_t z, uint64_t t) {
  const uint64_t degree = poly.get_length();
  const uint64_t mask = (1ULL << z) - 1;
  for (uint64_t i = 0; i < degree; ++i) {
    for (uint64_t j = 0; j < t; ++j) {
      mat[j][i] = (poly[i] >> (b + z * j)) & mask;
    }
  }
}

void ApproximateGadgetDecomp(RlweCiphertext& ct_bfv,
                             std::vector<Poly>& result_b,
                             std::vector<Poly>& result_a, uint64_t b,
                             uint64_t z, uint64_t t) {
  ApproximateGadgetDecomp(ct_bfv.b, result_b, b, z, t);
  ApproximateGadgetDecomp(ct_bfv.a, result_a, b, z, t);
}

}  // namespace

Poly::Poly() : length_(0), modulus_(0), nttform_(false) {}

Poly::Poly(uint64_t len, uint64_t modulus, bool ntt)
    : length_(len), modulus_(modulus), nttform_(ntt), payload(len, 0) {}

Poly::Poly(const std::vector<uint64_t>& vec, uint64_t modulus, bool ntt)
    : length_(vec.size()), modulus_(modulus), nttform_(ntt), payload(vec) {}

uint64_t Poly::get_length() const { return length_; }
uint64_t Poly::get_modulus() const { return modulus_; }
bool Poly::get_nttform() const { return nttform_; }
void Poly::set_nttform(bool ntt) { nttform_ = ntt; }
uint64_t* Poly::data() { return payload.data(); }
const uint64_t* Poly::data() const { return payload.data(); }
uint64_t& Poly::operator[](size_t idx) { return payload[idx]; }
const uint64_t& Poly::operator[](size_t idx) const { return payload[idx]; }
std::vector<uint64_t>& Poly::get_data() { return payload; }
const std::vector<uint64_t>& Poly::get_data() const { return payload; }

void PolyMult(Poly& lhs, Poly& rhs, Poly& result) {
  intel::hexl::EltwiseMultMod(result.data(), lhs.data(), rhs.data(),
                              result.get_length(), result.get_modulus(), 1);
}

void PolyAdd(Poly& lhs, Poly& rhs, Poly& result) {
  intel::hexl::EltwiseAddMod(result.data(), lhs.data(), rhs.data(),
                             result.get_length(), result.get_modulus());
}

void PolySub(Poly& lhs, Poly& rhs, Poly& result) {
  intel::hexl::EltwiseSubMod(result.data(), lhs.data(), rhs.data(),
                             result.get_length(), result.get_modulus());
}

void PolyFMA(Poly& mul, uint64_t constant, Poly& result) {
  intel::hexl::EltwiseFMAMod(result.data(), mul.data(), constant, nullptr,
                             result.get_length(), result.get_modulus(), 1);
}

void PolyToNTT(Poly& poly, intel::hexl::NTT& ntts) {
  if (poly.get_nttform()) {
    return;
  }
  ntts.ComputeForward(poly.data(), poly.data(), 1, 1);
  poly.set_nttform(true);
}

void PolyToCoef(Poly& poly, intel::hexl::NTT& ntts) {
  if (!poly.get_nttform()) {
    return;
  }
  ntts.ComputeInverse(poly.data(), poly.data(), 1, 1);
  poly.set_nttform(false);
}

Secret::Secret(std::vector<int64_t>& vec, uint64_t modulus, bool ntt) {
  std::vector<uint64_t> tmp(vec.size(), 0);
  for (size_t i = 0; i < vec.size(); ++i) {
    tmp[i] = (vec[i] < 0) ? (modulus - static_cast<uint64_t>(-vec[i]))
                          : static_cast<uint64_t>(vec[i]);
  }
  data_ = Poly(tmp, modulus, ntt);
  if (modulus == kLrCrtMod && vec.size() == 4096) {
    ntts_ = intel::hexl::NTT(vec.size(), modulus, kLrRootOfUnityCrt4096);
  } else if (modulus == kLrCrtMod && vec.size() == 2048) {
    ntts_ = intel::hexl::NTT(vec.size(), modulus, kLrRootOfUnityCrt2048);
  } else {
    ntts_ = intel::hexl::NTT(vec.size(), modulus);
  }
  if (ntt) {
    ntts_.ComputeForward(data_.data(), data_.data(), 1, 1);
    data_.set_nttform(true);
  }
}

intel::hexl::NTT& Secret::get_ntt() { return ntts_; }
Poly& Secret::get_data() { return data_; }
const Poly& Secret::get_data() const { return data_; }
uint64_t Secret::get_modulus() const { return data_.get_modulus(); }
uint64_t Secret::get_length() const { return data_.get_length(); }
bool Secret::get_nttform() const { return data_.get_nttform(); }
void Secret::to_ntt_form() { PolyToNTT(data_, ntts_); }

RlweCiphertext::RlweCiphertext() = default;

RlweCiphertext::RlweCiphertext(uint64_t len, uint64_t modulus, bool ntt)
    : a(len, modulus, ntt), b(len, modulus, ntt) {}

uint64_t RlweCiphertext::get_degree() const { return a.get_length(); }
uint64_t RlweCiphertext::get_modulus() const { return a.get_modulus(); }
void RlweCiphertext::set_nttform(bool ntt) {
  a.set_nttform(ntt);
  b.set_nttform(ntt);
}

Plaintext::Plaintext(uint64_t len, uint64_t pmod) : pt(len, pmod, false) {}

Plaintext::Plaintext(const std::vector<uint64_t>& vec, uint64_t pmod, bool ntt)
    : pt(vec, pmod, ntt) {}

uint64_t Plaintext::get_pmod() const { return pt.get_modulus(); }
bool Plaintext::get_isntt() const { return pt.get_nttform(); }

ApproximateRgswCiphertext::ApproximateRgswCiphertext(uint64_t len,
                                                     uint64_t modulus,
                                                     uint64_t b, uint64_t z,
                                                     uint64_t t)
    : degree_(len),
      modulus_(modulus),
      b_(b),
      z_(z),
      t_(t),
      ct_m(t, RlweCiphertext(len, modulus, true)),
      ct_sm(t, RlweCiphertext(len, modulus, true)) {}

uint64_t ApproximateRgswCiphertext::get_degree() const { return degree_; }
uint64_t ApproximateRgswCiphertext::get_modulus() const { return modulus_; }
uint64_t ApproximateRgswCiphertext::get_b() const { return b_; }
uint64_t ApproximateRgswCiphertext::get_z() const { return z_; }
uint64_t ApproximateRgswCiphertext::get_t() const { return t_; }

void RlweBfvEncrypt(Plaintext& pt, RlweCiphertext& ct, Secret& sk, float sig) {
  const uint64_t modulus = sk.get_modulus();
  const uint64_t length = sk.get_length();
  const uint64_t pmod = pt.get_pmod();
  const uint64_t delta = static_cast<uint64_t>(
      std::floor(static_cast<double>(modulus) / static_cast<double>(pmod)));

  Poly tmp(length, modulus, false);
  ct.set_nttform(true);
  SampleRandom(ct.a, modulus, length);

  if (!sk.get_nttform()) {
    sk.to_ntt_form();
  }

  intel::hexl::NTT& ntts = sk.get_ntt();
  if (pt.get_isntt()) {
    PolyToCoef(pt.pt, ntts);
  }

  PolyMult(ct.a, sk.get_data(), ct.b);
  PolyFMA(pt.pt, delta, tmp);

  Poly err(length, modulus, false);
  SampleGauss(err.get_data(), sig, modulus);

  PolyToCoef(ct.b, ntts);
  PolyAdd(ct.b, err, ct.b);
  PolyAdd(ct.b, tmp, ct.b);
  PolyToNTT(ct.b, ntts);
}

void RlweBfvEncode(Plaintext& pt, RlweCiphertext& ct, Secret& sk, float sig) {
  const uint64_t modulus = sk.get_modulus();
  const uint64_t length = sk.get_length();

  ct.set_nttform(true);
  SampleRandom(ct.a, modulus, length);
  if (!sk.get_nttform()) {
    sk.to_ntt_form();
  }

  PolyMult(ct.a, sk.get_data(), ct.b);
  intel::hexl::NTT& ntts = sk.get_ntt();
  if (pt.get_isntt()) {
    PolyToCoef(pt.pt, ntts);
  }

  Poly err(length, modulus, false);
  SampleGauss(err.get_data(), sig, modulus);

  PolyToCoef(ct.b, ntts);
  PolyAdd(ct.b, err, ct.b);
  PolyAdd(ct.b, pt.pt, ct.b);
  PolyToNTT(ct.b, ntts);
}

void RlweBfvDecrypt(RlweCiphertext& ct, Plaintext& pt, Secret& sk) {
  const uint64_t modulus = sk.get_modulus();
  const uint64_t pmod = pt.get_pmod();
  const uint64_t delta = static_cast<uint64_t>(
      std::floor(static_cast<double>(modulus) / static_cast<double>(pmod)));

  if (!sk.get_nttform()) {
    sk.to_ntt_form();
  }

  Poly tmp(ct.get_degree(), modulus, true);
  PolyMult(ct.a, sk.get_data(), tmp);
  PolySub(ct.b, tmp, tmp);
  intel::hexl::NTT& ntts = sk.get_ntt();
  PolyToCoef(tmp, ntts);

  for (uint64_t i = 0; i < ct.get_degree(); ++i) {
    int64_t centered = (tmp.payload[i] > (modulus >> 1))
                           ? static_cast<int64_t>(tmp.payload[i] - modulus)
                           : static_cast<int64_t>(tmp.payload[i]);
    int64_t decoded = static_cast<int64_t>(
        std::llround(static_cast<long double>(centered) /
                     static_cast<long double>(delta)));
    decoded %= static_cast<int64_t>(pmod);
    if (decoded < 0) {
      decoded += static_cast<int64_t>(pmod);
    }
    pt.pt[i] = static_cast<uint64_t>(decoded);
  }
}

void RlweBfvAdd(RlweCiphertext& add1, RlweCiphertext& add2,
                RlweCiphertext& result) {
  PolyAdd(add1.a, add2.a, result.a);
  PolyAdd(add1.b, add2.b, result.b);
}

void ExternalProduct(const ApproximateRgswCiphertext& ct_rgsw,
                     RlweCiphertext& ct_bfv, RlweCiphertext& result,
                     intel::hexl::NTT& ntts) {
  const uint64_t length = ct_bfv.get_degree();
  const uint64_t modulus = ct_bfv.get_modulus();
  const uint64_t t = ct_rgsw.get_t();

  std::vector<Poly> ct_bfv_b(t, Poly(length, modulus));
  std::vector<Poly> ct_bfv_a(t, Poly(length, modulus));

  ct_bfv.b.set_nttform(true);
  ct_bfv.a.set_nttform(true);
  PolyToCoef(ct_bfv.b, ntts);
  PolyToCoef(ct_bfv.a, ntts);
  ApproximateGadgetDecomp(ct_bfv, ct_bfv_b, ct_bfv_a, ct_rgsw.get_b(),
                          ct_rgsw.get_z(), t);

  for (uint64_t i = 0; i < t; ++i) {
    PolyToNTT(ct_bfv_b[i], ntts);
    PolyToNTT(ct_bfv_a[i], ntts);
  }

  RlweCiphertext tmp(length, modulus, true);
  RlweCiphertext result1(length, modulus, true);
  RlweCiphertext result2(length, modulus, true);
  for (uint64_t i = 0; i < t; ++i) {
    Poly ct_m_b = ct_rgsw.ct_m[i].b;
    Poly ct_m_a = ct_rgsw.ct_m[i].a;
    PolyMult(ct_bfv_b[i], ct_m_b, tmp.b);
    PolyMult(ct_bfv_b[i], ct_m_a, tmp.a);
    RlweBfvAdd(tmp, result1, result1);
  }
  for (uint64_t i = 0; i < t; ++i) {
    Poly ct_sm_b = ct_rgsw.ct_sm[i].b;
    Poly ct_sm_a = ct_rgsw.ct_sm[i].a;
    PolyMult(ct_bfv_a[i], ct_sm_b, tmp.b);
    PolyMult(ct_bfv_a[i], ct_sm_a, tmp.a);
    RlweBfvAdd(tmp, result2, result2);
  }
  PolySub(result1.a, result2.a, result.a);
  PolySub(result1.b, result2.b, result.b);
}

void RgswEncode(Plaintext& pt, ApproximateRgswCiphertext& ct, Secret sk,
                float sig) {
  const uint64_t z_gsw = ct.get_z();
  const uint64_t t_gsw = ct.get_t();
  const uint64_t b_gsw = ct.get_b();
  const uint64_t length = ct.get_degree();

  Poly sm(length, ct.get_modulus(), true);

  if (!sk.get_nttform()) {
    sk.to_ntt_form();
  }
  intel::hexl::NTT& ntts = sk.get_ntt();
  if (!pt.get_isntt()) {
    PolyToNTT(pt.pt, ntts);
  }
  PolyMult(pt.pt, sk.get_data(), sm);

  PolyToCoef(sm, ntts);
  PolyToCoef(pt.pt, ntts);

  RlweCiphertext tmp(length, ct.get_modulus(), true);
  PolyFMA(pt.pt, (1ULL << b_gsw), pt.pt);
  PolyFMA(sm, (1ULL << b_gsw), sm);
  const uint64_t base = (1ULL << z_gsw);
  for (uint64_t i = 0; i < t_gsw; ++i) {
    RlweBfvEncode(pt, tmp, sk, sig);
    tmp.b.set_nttform(true);
    tmp.a.set_nttform(true);
    ct.ct_m[i] = tmp;

    Plaintext pt_sm(sm.get_data(), ct.get_modulus(), false);
    RlweBfvEncode(pt_sm, tmp, sk, sig);
    tmp.b.set_nttform(true);
    tmp.a.set_nttform(true);
    ct.ct_sm[i] = tmp;

    PolyFMA(pt.pt, base, pt.pt);
    PolyFMA(sm, base, sm);
  }
}

}  // namespace psi::inspire::lrhe
