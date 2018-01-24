#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
namespace mp = boost::multiprecision;
using multiFloat = mp::number<mp::cpp_dec_float<400>>;

constexpr int dim = 2;
constexpr int planets = 3;
constexpr int N = planets * dim * 2;

//時刻に関するパラメータ
multiFloat dt("0.00000001");
//const multiFloat t_limit("1000000000000.0");
const multiFloat t_limit("100.10");

const multiFloat RTol("10e-12");
const multiFloat ATol("10e-12");
const multiFloat t_min("10e-100");
const multiFloat t_max("0.001");

//インターバル
constexpr int INTV = 0;
