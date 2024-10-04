
const Polynomial = require('polynomial');

Polynomial.setField({
    "add": function (a, b) {
        return BigInt(a) + BigInt(b);
    },
    "sub": function (a, b) {
        return BigInt(a) - BigInt(b);
    },
    "mul": function (a, b) {
        return BigInt(a) * BigInt(b);
    },
    "div": function (a, b) {
        return BigInt(a) / BigInt(b);
    },
    "parse": function (x) {
        return BigInt(x);
    },
    "empty": function (x) {
        return BigInt(x) === BigInt(0);
    },
    "pow": function (a, b) {
        return BigInt(a) ** BigInt(b);
    },
    "abs": function (a) {
        a = BigInt(a);
        if (a < 0n) {
            return -a;
        }
        return a;
    }
});

console.log(new Polynomial("3891238912389182391919282829x^3+4").add("23819218x^2-381").toString());