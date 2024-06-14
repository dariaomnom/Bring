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

DoubleDouble add(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble sum = utilTwoSum(a.hi, b.hi);
    double s = sum.hi;
    double e = sum.lo + a.lo + b.lo;
    sum = utilTwoSum(s, e);
    return sum;
}

DoubleDouble utilTwoProd(double a, double b) {
    double p = a * b;
    double err = std::fma(a, b, -p);
    return DoubleDouble(p, err);
}

DoubleDouble multiply(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble p1 = utilTwoProd(a.hi, b.hi);
    p1.lo += a.hi * b.lo + a.lo * b.hi;
    DoubleDouble prod = utilTwoSum(p1.hi, p1.lo);
    return prod;
}

DoubleDouble divide(const DoubleDouble& a, const DoubleDouble& b) {
    double q1 = a.hi / b.hi;
    DoubleDouble r = add(a, DoubleDouble(-q1 * b.hi, -q1 * b.lo));
    double q2 = r.hi / b.hi;
    DoubleDouble result = add(DoubleDouble(q1, 0), DoubleDouble(q2, 0));
    return result;
}

// f(x) = x^5 + x + a
DoubleDouble f(const DoubleDouble& x, double a) {
    DoubleDouble x2 = multiply(x, x);
    DoubleDouble x4 = multiply(x2, x2);
    DoubleDouble x5 = multiply(x4, x);
    DoubleDouble result = add(x5, DoubleDouble(x.hi, x.lo));
    return add(result, DoubleDouble(a, 0));
}

// Производная f(x)
DoubleDouble f_prime(const DoubleDouble& x) {
    DoubleDouble x2 = multiply(x, x);
    DoubleDouble x4 = multiply(x2, x2);
    DoubleDouble five_x4 = multiply(DoubleDouble(5.0, 0.0), x4);
    return add(five_x4, DoubleDouble(1.0, 0.0));
}

// Метод Ньютона
DoubleDouble newton_method(double a, DoubleDouble x0, int max_iter = 1000, double tol = 1e-32) {
    for (int i = 0; i < max_iter; ++i) {
        DoubleDouble fx = f(x0, a);
        DoubleDouble fx_prime = f_prime(x0);

        DoubleDouble dx = divide(DoubleDouble(-fx.hi, -fx.lo), fx_prime);

        x0 = add(x0, dx);

        printf("Iter %d: \n x = %.40e   %.40e\n", i+1, x0.hi, x0.lo);
        printf(" F(x) = %.40e   %.40e \n", fx.hi, fx.lo);
        printf(" dx/x = %.40e   %.40e \n", dx.hi / x0.hi, dx.lo / x0.lo);

        // Условие остановки на основе относительного изменения
        if (std::fabs(dx.hi / x0.hi) < tol && std::fabs(dx.lo / x0.lo) < tol) {
            break;
        }
    }
    return x0;
}

// Метод хорд
DoubleDouble secant_method(double a, DoubleDouble x0, DoubleDouble x1, int max_iter = 1000, double tol = 1e-32) {
    for (int i = 0; i < max_iter; ++i) {
        DoubleDouble fx0 = f(x0, a);
        DoubleDouble fx1 = f(x1, a);

        DoubleDouble denom = add(fx1, DoubleDouble(-fx0.hi, -fx0.lo));
        DoubleDouble num = add(x1, DoubleDouble(-x0.hi, -x0.lo));
        DoubleDouble slope = divide(denom, num);
        DoubleDouble dx = divide(DoubleDouble(-fx1.hi, -fx1.lo), slope);

        x0 = x1;
        x1 = add(x1, dx);

        printf("Iter %d: \n x = %.40e   %.40e\n", i+1, x1.hi, x1.lo);
        printf(" from F(x) = %.40e   %.40e \n to   F(x) = %.40e   %.40e\n", fx0.hi, fx0.lo, fx1.hi, fx1.lo);
        printf(" dx/x = %.40e   %.40e\n", dx.hi / x1.hi, dx.lo / x1.lo);

        if (std::fabs(dx.hi / x1.hi) < tol && std::fabs(dx.lo / x1.lo) < tol) {
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
