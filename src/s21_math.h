#ifndef S21_MATH_H_
#define S21_MATH_H_

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define s21_PI 3.141592653589793
#define s21_EPS 1E-16
#define s21_INF __builtin_inf()
#define s21_NAN __builtin_nan("")
#define s21_NEG_ZERO -0.0
#define ACCURACY 100

#define TEST_EPS 1E-6  //для тестов

/*----Модули----*/
int s21_abs(int x);
long double s21_fabs(double x);

/*----Тригонометрические функции----*/
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_cos(double x);
long double s21_sin(double x);
long double s21_tan(double x);

/*----Степень, экспонента, логарифм----*/
long double s21_exp(double x);  // экспонента в степени х
long double s21_log(double x);  // натуральный логарифм
long double s21_pow(double base, double exp);
long double s21_sqrt(double x);

/*----Целочсленное округление + остаток от деления----*/
long double s21_ceil(double x);   // округление вверх
long double s21_floor(double x);  // округление вниз
long double s21_fmod(double x,
                     double y);  // остаток от деления с плавающей точкой

/*----Вспомогательные функции----*/
int s21_isnan(long double x);
int s21_isinf(long double x);
long double s21_factorial(double x);

#endif