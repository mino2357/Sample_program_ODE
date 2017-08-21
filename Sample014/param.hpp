#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
namespace mp = boost::multiprecision;
using cpp_dec_float = mp::number<mp::cpp_dec_float<50>>;
using multiFloat = cpp_dec_float;

constexpr int dim = 2;
constexpr int planets = 10;
constexpr int N = planets * dim * 2;

//時刻に関するパラメータ
multiFloat dt("0.0001");
const multiFloat t_limit("1000000000000.0");

const multiFloat RTol("10e-10");
const multiFloat ATol("10e-10");
const multiFloat t_min("10e-50");
const multiFloat t_max("0.01");

//インターバル
constexpr int INTV = 1;
