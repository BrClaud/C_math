#include "s21_math.h"

/*----Модуль числа с плавающей точкой----*/
long double s21_fabs(double x) {
  if (x < 0.0) x = -x;
  return (long double)x;
}

/*----Модуль целого числа----*/
int s21_abs(int x) {
  if (x < 0) x = -x;
  return (int)x;
}

/*----Возвращает 1 если число nan----*/
int s21_isnan(long double x) { return (x != x) ? 1 : 0; }

/*----Возвращает 1 если число бесконечность----*/
int s21_isinf(long double x) {
  return ((x == s21_INF || x == -s21_INF) && !s21_isnan(x)) ? 1 : 0;
}

/*----Возвращает ближайшее целое число, не превышающее заданное значение----*/
long double s21_floor(double x) {
  // округляем если больше 0
  long double res = (long long int)x;
  if (s21_isinf(x) || x == 0.0 || x == s21_NEG_ZERO || s21_isnan(x) ||
      res == LONG_MIN || res == LONG_MAX || x == DBL_MAX) {
    res = x;
  } else {
    // округляем если меньше 0
    if (x < 0.0 && x != (double)res) res--;
  }
  return res;
}

/*----Возвращает ближайшее целое число, не меньшее заданного значения----*/
long double s21_ceil(double x) {
  long double res = (long long int)x;
  if (s21_isinf(x) || x == s21_NEG_ZERO || s21_isnan(x) || res == LONG_MIN ||
      res == LONG_MAX || x == DBL_MAX) {
    res = x;
  } else {
    res = s21_floor(x) + (x != (long long int)x);
  }
  return res;
}

/*----Остаток от деления x на y----*/
/* Функция рассчитывает остаток от деления x на y по формуле x-n*y,
где n – это результат деления x на y округленный в сторону нуля */
long double s21_fmod(double x, double y) {
  long double res = 0.0;
  // если х или у нан то возвращаем нан, если х бесконечность или если делим на
  // 0 то возвращаем нан
  if (s21_isnan(x) || s21_isnan(y) || s21_isinf(x) || y == 0) {
    res = s21_NAN;
  }
  //если делим на бесконечность то возращаем х без изменений
  else if (s21_isinf(y)) {
    res = x;
  } else {
    res = (int)(x / y);
    res = x - res * y;
  }
  return res;
}

// add_value это то что прибавляем (новый член ряда)
// res это наша общая сумма ряда которую мы считаем
/*e^x = 1 + x + x^2/2 + x^3/6 + x^4/24 + ...
e^x = 1 + x + (x^2)/(2!) + (x^3)/(3!) + (x^4)/(4!) + ...*/
long double s21_exp(double x) {
  long double add_value = 1.0, res = 1.0, i = 1.0;
  if (x == s21_INF) {
    res = s21_INF;
  } else if (s21_isnan(x)) {
    res = s21_NAN;
  } else if (x == -s21_INF) {
    res = 0.0;
  } else {
    long double copy_x = x;
    if (x < 0.0) copy_x = -x;
    int flag = 1;
    // так как каждый слудеющий член ряда мы делим на (i+1)! то он уменьшается с
    // каждый итерацией и однажды он точно будет меньше эпсилона
    while (s21_fabs(add_value) > s21_EPS && flag) {
      add_value *=
          copy_x / i;  // увеличиваем степень икса и домножаем факториал
                       // (было 1*2*3 = 6 это 3!, стало 6*4 = 24 это 4!)
      i++;
      res += add_value;
      if (res > DBL_MAX) {
        res = s21_INF;
        flag = 0;
      }
    }
    if (x < 0.0) res = 1 / res;
    /* Когда x отрицательное, мы вычисляем модуль, а экспонента в степени
    отрицательного числа будет равна обратному значению экспоненты от его
    модуля */
  }
  return res;
}

/*----Вычисляет натуральный логарифм----*/
// метод Галея
long double s21_log(double x) {
  long double res = 0;
  if (x == s21_INF) {
    res = s21_INF;
  } else if (x == 0) {
    res = -s21_INF;
  } else if (x < 0) {
    res = s21_NAN;
  } else {
    long double previous = 0;
    int e_repeat = 0;
    // e_repeat используется для подсчёта того, сколько раз нужно нужно
    // разделить x, чтобы получить значение, меньшее, чем экспонента. Это нужно
    // для обеспечения лучшей точности (как в тригонометрии загоняются значения
    // от -2пи до 2 пи так и тут в промежуток меньший чем экспонента).
    for (; x >= s21_exp(1); e_repeat++) x /= s21_exp(1);
    for (int i = 0; i < ACCURACY; i++) {
      previous = res;
      res = previous + 2 * (x - s21_exp(previous)) / (x + s21_exp(previous));
    }
    // корректируем результат на количество степеней, на которое сократили число
    res += e_repeat;
  }
  return res;
}

/*----Возводит число в заданную степень----*/
// Формула: b^a = e^(a*ln(b))
/* Для отрицательных числел степень должна быть целая, т.к.
 (a^(p/q) где p - целое, а q - натуральное, определена только
 при a>0, т.к. (-8)^(1/3)!=(-8)^(2/6) - возникает неопределённость, в
 первом случае -2, а в 2, xотя 1/3 = 2/6) */
long double s21_pow(double base, double exp) {
  long double res = 0;
  if (exp == 0)  //возведение в нулевую степень
    res = 1;
  else if (s21_isnan(base) || s21_isnan(exp))
    res = s21_NAN;
  else if (base == 0) {  //возведение нуля в степень
    if (exp > 0)
      res = base;
    else
      res = s21_INF;
  } else if (s21_isinf(base) && s21_isinf(exp))
    res = (exp < 0) ? 0 : s21_INF;
  else if ((base == -1 || base == 1) && s21_isinf(exp))
    res = 1.0;  // краевые случаи
  else {
    if (base > 0)  //основание положительное
      res = s21_exp(exp * s21_log(base));
    else if (s21_ceil(exp) ==
             exp)  // основание отрицательное, надо чтобы степень была целая
    {
      if (s21_abs(((long long int)exp) % 2) ==
          1)  // степень - нечётная => работаем с модулем и возвращаем минус
        res = -s21_exp(exp * s21_log(s21_fabs(base)));
      else  // степень - чётная => просто работаем с модулем
        res = s21_exp(exp * s21_log(s21_fabs(base)));
    } else if (base == -s21_INF) {  // основание отрицательное, но у минус
                                    // бесконечности свои правила
      if (exp < 0) res = 0.0;
      if (exp > 0) res = s21_INF;
    } else  // в ином случае нан
      res = s21_NAN;
  }
  return res;
}

/*----Вычисляет квадратный корень----*/
long double s21_sqrt(double x) { return s21_pow(x, 0.5); }

/*----Вычисляет синус----*/
// ряд Тейлора (аппроксимация значения синуса в окрестности нуля)
long double s21_sin(double x) {
  long double res = 0.0;
  double add_value = 1.0;
  int i = 0;
  if (s21_isinf(x) || s21_isnan(x))
    res = s21_NAN;
  else {
    // приводим x к диапазону до 2pi
    x = s21_fmod(x, 2 * s21_PI);
    // формула ряда
    while (s21_fabs(add_value) >= s21_EPS) {
      add_value =
          s21_pow(-1, i) * s21_pow(x, 2 * i + 1) / s21_factorial(2 * i + 1);
      res += add_value;
      i++;
    }
  }
  return res;
}

/*----Вычисляет факториал----*/
// 4! = 1*2*3*4
long double s21_factorial(double x) {
  long double res = 1.0;
  if (x >= 0)
    for (int i = 1; i <= x; i++) res *= i;
  else
    res = 0;
  return res;
}

/*----Вычисляет косинус----*/
long double s21_cos(double x) { return s21_sin(s21_PI / 2 + x); }

/*----Вычисляет тангенс----*/
long double s21_tan(double x) {
  long double res = s21_sin(x) / s21_cos(x);
  return res;
}

/*----Вычисляет арктангенс----*/
// ряд Маклорена  arctan(x) = x - (x^3)/3 + (x^5)/5 - (x^7)/7 + ...
long double s21_atan(double x) {
  long double res = x, add_value = x, i = 1.0;
  if (x == 1) {
    res = s21_PI / 4;
  } else if (x == -1) {
    res = -s21_PI / 4;
  } else if ((x < 1) && (x > -1)) {
    while (s21_fabs(add_value) > s21_EPS) {
      // -1.0 используется для чередования знаков при сложении и вычитании
      // членов ряда Маклорена
      add_value = -1.0 * add_value * x * x * (2.0 * i - 1.0) / (2.0 * i + 1.0);
      res += add_value;
      i++;
    }
  } else if (x > 1) {  // 1>0
    res = s21_PI / 2 - s21_atan(1 / x);
  } else if (x < -1) {  //-1<0
    res = -s21_PI / 2 - s21_atan(1 / x);
  }
  return res;
}

/*----Вычисляет арксинус----*/
long double s21_asin(double x) {
  long double res = 0;
  if (x < -1 || x > 1 || s21_isnan(x))
    res = s21_NAN;
  else if (x == 1)
    res = s21_PI / 2;
  else if (x == -1)
    res = -s21_PI / 2;
  else
    res = s21_atan(x / s21_sqrt(1 - x * x));
  return res;
}

/*----Вычисляет арккосинус----*/
long double s21_acos(double x) { return s21_PI / 2 - s21_asin(x); }