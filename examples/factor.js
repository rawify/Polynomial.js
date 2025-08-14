
var Polynomial = require('polynomial');

Polynomial.setField("Q");

// This implements Newton's method on top of a rational polynomial. The code is highly experimental!
function factor(P) {

    var xn = 1;

    var F = Polynomial(P); // f(x)
    var f = F.derive(); // f'(x)

    var res = [];

    do {

        var prev, xn = Polynomial(1).coeff[0];

        do {
            prev = xn;
            xn = xn.sub(F.eval(xn).div(f.eval(xn)));
        } while (Math.abs(xn.valueOf() - prev.valueOf()) > 1e-10);

        var p = Polynomial("x").sub(xn); // x-x0

        F = F.div(p);
        f = F.derive();             

        res.push(xn);

    } while (f.degree() >= 0);

    return res;
}
// Result should by ~3, ~4, ~9
console.log(factor("x^3-16x^2+75x-108").map(x => x.toString()));