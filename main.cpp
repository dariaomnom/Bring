#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

struct DoubleDouble {
    double hi;
    double lo;
    explicit DoubleDouble(double h = 0, double l = 0) : hi(h), lo(l) {}
};

DoubleDouble utilTwoSum(double a, double b) {
    double s = a + b;
    double t = s - a;
    double err = (a - (s - t)) + (b - t);
    return DoubleDouble(s, err);
}

//DoubleDouble add(const DoubleDouble& a, const DoubleDouble& b) {
//    DoubleDouble s1 = utilTwoSum(a.hi, b.hi);
//    DoubleDouble s2 = utilTwoSum(a.lo, b.lo);
//    return utilTwoSum(s1.hi, s2.hi);
//}

//DoubleDouble add(const DoubleDouble& a, const DoubleDouble& b) {
//    DoubleDouble s1 = utilTwoSum(a.hi, b.hi);
//    DoubleDouble s2 = utilTwoSum(a.lo, b.lo);
//    return utilTwoSum(s1.hi, s1.lo + s2.hi);
//}

DoubleDouble utilThreeSum(double a, double b, double c) {
    DoubleDouble sum1 = utilTwoSum(a, b);
    DoubleDouble sum2 = utilTwoSum(sum1.hi, c);
    double err = sum1.lo + sum2.lo;
    return DoubleDouble(sum2.hi, err);
}

DoubleDouble add(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble s1 = utilTwoSum(a.hi, b.hi);
    DoubleDouble s2 = utilTwoSum(a.lo, b.lo);
    return utilThreeSum(s1.hi, s1.lo, s2.hi);
}

DoubleDouble utilTwoProd(double a, double b) {
    double p = a * b;
    double err = std::fma(a, b, -p);
    return DoubleDouble(p, err);
}

DoubleDouble multiply(const DoubleDouble& a, const DoubleDouble& b) {
//    DoubleDouble p1 = utilTwoProd(a.hi, b.hi);
//    DoubleDouble p2 = utilTwoSum(a.hi * b.lo, a.lo * b.hi);
//    DoubleDouble p3 = utilTwoSum(p1.lo, p2.hi);
//    return utilTwoSum(p1.hi, p3.hi);
    DoubleDouble p1 = utilTwoProd(a.hi, b.hi);
    DoubleDouble p2 = utilTwoProd(a.hi, b.lo);
    DoubleDouble p3 = utilTwoProd(a.lo, b.hi);
    DoubleDouble sum1 = add(p1, DoubleDouble(p2.hi + p3.hi, p2.lo + p3.lo));
    return add(sum1, DoubleDouble(0, p1.lo));
}

//DoubleDouble divide2(const DoubleDouble& a, const DoubleDouble& b) {
//    double q1 = a.hi / b.hi;
//    DoubleDouble prod1 = multiply(DoubleDouble(q1, 0.0), b);
//    DoubleDouble r = add(a, DoubleDouble(-prod1.hi, -prod1.lo));
//    double q2 = r.hi / b.hi;
//    DoubleDouble prod2 = multiply(DoubleDouble(q2, 0.0), b);
//    r = add(r, DoubleDouble(-prod2.hi, -prod2.lo));
//    double q3 = r.hi / b.hi;
//    return add(DoubleDouble(q1, 0.0), add(DoubleDouble(q2, 0.0), DoubleDouble(q3, 0.0)));
//}

DoubleDouble divide(const DoubleDouble& a, const DoubleDouble& b) {
    double q1 = a.hi / b.hi;
    DoubleDouble r = add(a, DoubleDouble(-q1 * b.hi, -q1 * b.lo));
    double q2 = r.hi / b.hi;
    DoubleDouble result = add(DoubleDouble(q1, 0), DoubleDouble(q2, 0));
    return result;
}

static DoubleDouble subtract(const DoubleDouble& x, const DoubleDouble& y){
    DoubleDouble S = utilTwoSum(x.hi, -y.hi);
    DoubleDouble E = utilTwoSum(x.lo, -y.lo);
    double c = S.lo + E.hi;
    double vh = S.hi + E.hi;
    double vl = c - (vh - S.hi);
    c = vl + E.lo;
    return DoubleDouble(vh + c, c - (vh + c - vh));
}

//static DoubleDouble divide(const DoubleDouble& a, const DoubleDouble& b) {
//    double q1 = a.hi / b.hi;
//    DoubleDouble p = multiply(DoubleDouble(q1), b);
//    DoubleDouble r = subtract(a, DoubleDouble(p.hi, p.lo));
//    double q2 = r.hi / b.hi;
//    return add(DoubleDouble(q1), DoubleDouble(q2));
//}

// f(x) = x^5 + x + a
DoubleDouble f(const DoubleDouble& x, double a) {
    DoubleDouble x2 = multiply(x, x);
    DoubleDouble x4 = multiply(x2, x2);
    DoubleDouble x5 = multiply(x4, x);
    DoubleDouble ax = utilTwoSum(a, x.hi);
    return add(x5, DoubleDouble(ax.hi, ax.lo + x.lo));
}

// Производная f(x)
DoubleDouble f_prime(const DoubleDouble& x) {
    DoubleDouble x2 = multiply(x, x);
    DoubleDouble x4 = multiply(x2, x2);
    DoubleDouble five_x4 = multiply(DoubleDouble(5.0, 0.0), x4);
    return add(five_x4, DoubleDouble(1.0, 0.0));
}

// Метод Ньютона
DoubleDouble newton_method(double a, DoubleDouble x0, int max_iter = 1000, double tol = 1e-20) {
    for (int i = 0; i < max_iter; ++i) {
        DoubleDouble fx = f(x0, a);
        DoubleDouble fx_prime = f_prime(x0);

//        double dx_hi = -fx.hi / fx_prime.hi;
//        double dx_lo = (-fx.lo / fx_prime.hi) - (fx.hi * fx_prime.lo / (fx_prime.hi * fx_prime.hi));
//        DoubleDouble dx(dx_hi, dx_lo);
        DoubleDouble dx = divide(DoubleDouble(-fx.hi, -fx.lo), fx_prime);

        x0 = add(x0, dx);
        printf("Iter %d: \n x = %.16e   %.16e\n", i+1, x0.hi, x0.lo);
        printf(" F(x) = %.16e   %.16e \n", fx.hi, fx.lo);
//        printf("Newton Iteration %d: x = %.40e   %.40e\n", i+1, x0.hi, x0.lo);
        if (std::fabs(dx.hi) < tol && std::fabs(dx.lo) < tol) {
            break;
        }
    }
    return x0;
}

// Метод хорд
DoubleDouble secant_method(double a, DoubleDouble x0, DoubleDouble x1, int max_iter = 1000, double tol = 1e-20) {
    for (int i = 0; i < max_iter; ++i) {
        DoubleDouble fx0 = f(x0, a);
        DoubleDouble fx1 = f(x1, a);
//        printf("Iter %d: \n from F(x) = %.40e   %.40e \n to   F(x) = %.40e   %.40e\n",  i+1, fx0.hi, fx0.lo, fx1.hi, fx1.lo);

        // denominator = f(x1) - f(x0)
//        DoubleDouble denom = add(fx1, DoubleDouble(-fx0.hi, -fx0.lo));

        // numerator = x1 - x0
//        DoubleDouble num = add(x1, DoubleDouble(-x0.hi, -x0.lo));

        // slope = denom / num
//        DoubleDouble slope = divide(denom, num);
//        double slope_hi = denom.hi / num.hi;
//        double slope_lo = (denom.lo / num.hi) - (denom.hi * num.lo / (num.hi * num.hi));

        // dx = f(x1) / slope
//        DoubleDouble dx = divide(DoubleDouble(-fx1.hi, -fx1.lo), slope);
//        double dx_hi = -fx1.hi / slope_hi;
//        double dx_lo = (-fx1.lo / slope_hi) - (fx1.hi * slope_lo / (slope_hi * slope_hi));
//        DoubleDouble dx(dx_hi, dx_lo);

        DoubleDouble denom = add(fx1, DoubleDouble(-fx0.hi, -fx0.lo));
        DoubleDouble num = add(x1, DoubleDouble(-x0.hi, -x0.lo));
        DoubleDouble slope = divide(denom, num);
        DoubleDouble dx = divide(DoubleDouble(-fx1.hi, -fx1.lo), slope);

        x0 = x1;
        x1 = add(x1, dx);

//        printf("Secant Iteration %d: x = %.40e   %.40e\n", i+1, x1.hi, x1.lo);
        printf("Iter %d: \n x = %.16e   %.16e\n", i+1, x1.hi, x1.lo);
        printf(" from F(x) = %.16e   %.16e \n to   F(x) = %.16e   %.16e\n", fx0.hi, fx0.lo, fx1.hi, fx1.lo);
//        printf("Iter %d: \n from F(x) = %.40e   %.40e \n to   F(x) = %.40e   %.40e\n",  i+1, fx0.hi, fx0.lo, fx1.hi, fx1.lo);
        if (std::fabs(dx.hi) < tol && std::fabs(dx.lo) < tol) {
            break;
        }
    }
    return x1;
}

int main() {
    double a;
    std::cout << "x^5 + x + a = 0" << std::endl;
    std::cout << "Enter the value of a: ";
    std::cin >> a;

    double left, right;

    left = -1 * copysign(1, a) * pow(abs(a), 1.0/5) - 1;
    right = -1 * copysign(1, a) * pow(abs(a), 1.0/5) + 1;

    DoubleDouble x0(left, 0.0); // Начальное приближение для метода Ньютона
    DoubleDouble x1(right, 0.0); // Начальное приближение для метода хорд

    printf("\nInitial approximations: from %.1f to %.1f\n", left, right);


    std::cout << "\nUsing Newton's method:\n";
    DoubleDouble root_newton = newton_method(a, x0);

    std::cout << "\nUsing Secant method:\n";
    DoubleDouble root_secant = secant_method(a, x0, x1);

    printf("\nRoot found by Newton's method:  %.40e   %.40e\n", root_newton.hi, root_newton.lo);
    printf("Root found by Secant method:    %.40e   %.40e\n", root_secant.hi, root_secant.lo);

    return 0;
}
