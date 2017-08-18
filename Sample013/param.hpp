#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
namespace mp = boost::multiprecision;
using multiFloat = mp::cpp_dec_float_100;


constexpr int dim = 2;
constexpr int planets = 3;
constexpr int N = planets * dim * 2;

//時刻に関するパラメータ
multiFloat dt("1.0e-1");
const multiFloat t_limit("200.0");

const multiFloat RTol("10e-20");
const multiFloat ATol("10e-20");
const multiFloat t_min("10e-50");
const multiFloat t_max("0.1");

//インターバル
constexpr int INTV = 1;
