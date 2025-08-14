/**
 * @license Polynomial.js v1.6.0 8/14/2025
 * https://github.com/rawify/Polynomial.js
 *
 * Copyright (c) 2025, Robert Eisele (https://raw.org/)
 * Licensed under the MIT license.
 **/

import Fraction from 'fraction.js';
import Complex from 'complex.js';
import Quaternion from 'quaternion';


const STR_REGEXP = /([+-]?)(?:([^+x-]+)?(?:x(?:\^([\d\/]+))?)|([^+x-]+))/g;

/**
 * The actual field selected
 * 
 * @type Object
 */
const FIELDS = {
  'R': { // Run in R
    'add': function (a, b) { return a + b; },
    'sub': function (a, b) { return a - b; },
    'mul': function (a, b) { return a * b; },
    'div': function (a, b) {
      if (b === 0) throw new Error('DIV/0');
      return a / b;
    },
    'parse': function (x) { return parseFloat(x); },
    'empty': function (x) { return x === undefined || x === 0; },
    'pow': function (a, b) { return Math.pow(a, b); },
    'equals': function (a, b) { return a === b; }
  }
};

let FIELD = FIELDS['R'];


/**
 * Polynomial constructor
 * @constructor
 * @param {String|Object|number|null} x
 */
function Polynomial(x) {

  if (!(this instanceof Polynomial)) {
    return new Polynomial(x);
  }
  this['coeff'] = parse(x);
}

// Trace poly div steps
Polynomial['trace'] = null;

/**
 * Modular inverse for integers
 * @param {number} z
 * @param {number} n
 * @returns {number}
 */
const modinv = function (z, n) {

  /**
   *    z * s + n * t = 1 
   * => z * s mod n = 1
   * => z^-1 = s mod n
   */
  const tmp = egcd(z, n);
  if (tmp[0] !== 1) {
    throw new Error('DIV/-');
  }
  return mod(tmp[1], n);
};

/**
 * Calculates the gcd
 * 
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
function gcd(a, b) {

  for (let t; b;) {
    t = a;
    a = b;
    b = t % b;
  }
  return Math.abs(a);
}

function clone(x) {
  return Object.assign({}, x);
}

/**
 * Calculates the extended gcd
 * 
 * @param {number} a
 * @param {number} b
 * @returns {Array}
 */
const egcd = function (a, b) {

  // gcd = a * s  +  b * t

  let s = 0, t = 1, u = 1, v = 0;
  while (a !== 0) {

    const q = b / a | 0, r = b % a;
    const m = s - u * q, n = t - v * q;

    b = a;
    a = r;
    s = u;
    t = v;
    u = m;
    v = n;
  }
  return [b /* gcd*/, s, t];
};

/**
 * Mathematical modulo (always non-negative)
 * @param {number} n
 * @param {number} m
 * @returns {number}
 */
const mod = function (n, m) {
  return (n % m + m) % m;
};

/**
 * Falling factorial (n)_k = n * (n-1) * ... * (n-k+1)
 * with (n)_0 = 1
 * @param {number} n
 * @param {number} k
 * @returns {number}
 */
const fallingFactorial = function (n, k) {
  let p = 1;
  for (; k-- > 0; p *= n--){}
  return p;
};

/**
 * Rising factorial n^(k) = (n+1)(n+2)...(n+k)
 * with n^(0) = 1
 * @param {number} n
 * @param {number} k
 * @returns {number}
 */
const risingFactorial = function (n, k) {
  let p = 1;
  for (; k-- > 0; p *= ++n){}
  return p;
};

/**
 * Combines the keys of two objects
 * 
 * @param {Object} a
 * @param {Object} b
 * @returns {Object}
 */
function keyUnion(a, b) {

  const k = {};
  for (let i in a) {
    k[i] = 1;
  }
  for (let i in b) {
    k[i] = 1;
  }
  return k;
}

/**
 * Ensure we only keep non-empty coefficients; coerce keys to numbers.
 * Mutates and returns the same object for chaining.
 * @param {Object} obj
 * @returns {Object}
 */
const normalize = function (obj) {
  for (let k in obj) {
    if (FIELD['empty'](obj[k])) delete obj[k];
  }
  return obj;
};

/**
 * Gets the degree of the actual polynomial
 * 
 * @param {Object} x
 * @returns {number}
 */
function degree(x) {

  let i = -Infinity;

  for (let k in x) {
    if (!FIELD['empty'](x[k]))
      i = Math.max(+k, i);
  }
  return i;
}

/**
 * Division helper: given numerator x (mutated) and denominator y, return quotient q.
 * Also records steps if Polynomial['trace'] is set.
 * 
 * @param {Object} x The numerator coefficients
 * @param {Object} y The denominator coefficients
 * @returns {Object}
 */
const div = function (x, y) {

  const r = {};

  let i = degree(x);
  const j = degree(y);
  const trace = [];

  if (j < 0) {
    throw new Error('DIV/0');
  }

  while (i >= j) {

    const tmp = r[i - j] = FIELD['div'](x[i] || 0, y[j] || 0);

    for (let k in y) {
      x[+k + i - j] = FIELD['sub'](x[+k + i - j] || 0, FIELD['mul'](y[k] || 0, tmp));
    }

    if (Polynomial['trace'] !== null) {

      const tr = {};
      for (let   k in y) {
        tr[+k + i - j] = FIELD['mul'](y[k] || 0, tmp);
      }
      trace.push(new Polynomial(tr));
    }

    i = degree(x);
  }

  // Add rest
  if (Polynomial['trace'] !== null) {
    trace.push(new Polynomial(x));
    Polynomial['trace'] = trace;
  }
  return normalize(r);
};

/**
 * Leading exponent key (numeric) or null if empty
 * @param {Object} poly
 * @returns {number|null}
 */
function lcExp(poly) {

  let max = null;

  for (let i in poly) {

    if (!FIELD['empty'](poly[i])) {

      if (max === null || +max < +i) {
        max = i;
      }
    }
  }
  return +max;
}

/**
 * Make polynomial monic: divide by leading coefficient (in-place)
 * @param {Object} a
 * @param {number|null} max leading exponent (if known)
 * @returns {Object}
 */
function makeMonic(a, max) {

  if (max !== null) {

    for (let i in a) {
      a[i] = FIELD['div'](a[i], a[max]);
    }
  }
  return a;
}

/**
 * Parse multiplicative expression part (handles a*b/c chains)
 */
function parseExp(sgn, exp) {

  exp = String(exp).match(/[^*/]+|[*/]/g);

  let num = FIELD['parse'](sgn + exp[0]);

  for (let i = 1; i < exp.length; i += 2) {

    if (exp[i] === '*') {
      num = FIELD['mul'](num, FIELD['parse'](exp[i + 1] || 1));
    } else if (exp[i] === '/') {
      num = FIELD['div'](num, FIELD['parse'](exp[i + 1] || 1));
    }
  }
  return num;
}

/**
 * Parse polynomial input
 * 
 * @param {String|Object|number|null} x
 * @returns {Object} sparse coefficient map
 */
const parse = function (x) {

  const ret = {};

  if (x === null || x === undefined) {
    x = 0;
  }

  switch (typeof x) {

    case "object":

      if (x['coeff']) {
        x = x['coeff'];
      }

      if (Fraction && x instanceof Fraction || Complex && x instanceof Complex || Quaternion && x instanceof Quaternion) {
        ret[0] = x;
      } else
        // Handles Arrays the same way
        for (let i in x) {

          if (!FIELD['empty'](x[i])) {
            ret[i] = FIELD['parse'](x[i]);
          }
        }
      return ret;

    case "number":
      ret[0] = FIELD['parse'](x);
      return ret;

    case "string":

      const s = x.replace(/\s+/g, '');
      STR_REGEXP.lastIndex = 0;

      for (let m; (m = STR_REGEXP['exec'](s)) !== null;) {
        let num = 1, exp = 1;

        if (m[4] !== undefined) {
          num = m[4];
          exp = 0;
        } else if (m[2] !== undefined) {
          num = m[2];
        }

        num = parseExp(m[1], num);

        // Parse exponent
        if (m[3] !== undefined)
          exp = parseInt(m[3], 10);

        ret[exp] = ret[exp] === undefined
          ? num
          : FIELD['add'](ret[exp], num);
      }
      return ret;
  }
  throw new Error("Invalid Param");
};

/**
 * Helper method to stringify
 * 
 * @param {string} fn the callback name
 * @returns {Function}
 */
const toString = function (fn) {

  /**
   * The actual to string function 
   * 
   * @returns {string|null}
   */
  const Str = function () {

    const poly = this['coeff'];

    const keys = [];
    for (const i in poly) {
      keys.push(+i);
    }

    if (keys.length === 0)
      return "0";

    keys.sort(function (a, b) {
      return a - b;
    });

    let str = "";
    for (let k = keys.length; k--;) {

      const i = keys[k];

      const cur = poly[i];

      let val = cur;

      if (val === null || val === undefined)
        continue;

      if (Complex && val instanceof Complex) {

        // Add real part
        if (val['re'] !== 0) {

          if (str !== "" && val['re'] > 0) {
            str += "+";
          }

          if (val['re'] === -1 && i !== 0) {
            str += "-";
          } else if (val['re'] !== 1 || i === 0) {
            str += val['re'];
          }

          // Add exponent if necessary, no DRY, let's feed gzip
          if (i === 1)
            str += "x";
          else if (i !== 0)
            str += "x^" + i;
        }

        // Add imaginary part
        if (val['im'] !== 0) {

          if (str !== "" && val['im'] > 0) {
            str += "+";
          }

          if (val['im'] === -1) {
            str += "-";
          } else if (val['im'] !== 1) {
            str += val['im'];
          }

          str += "i";

          // Add exponent if necessary, no DRY, let's feed gzip
          if (i === 1)
            str += "x";
          else if (i !== 0)
            str += "x^" + i;
        }

      } else {

        val = val.valueOf();

        // Skip if it's zero
        if (val === 0)
          continue;

        // Separate by +
        if (str !== "" && val > 0) {
          str += "+";
        }

        if (val === -1 && i !== 0)
          str += "-";
        else

          // Add number if it's not a "1" or the first position
          if (val !== 1 || i === 0)
            str += cur[fn] ? cur[fn]() : cur['toString']();

        // Add exponent if necessary, no DRY, let's feed gzip
        if (i === 1)
          str += "x";
        else if (i !== 0)
          str += "x^" + i;
      }
    }

    if (str === "")
      return "0";

    return str;
  };
  return Str;
};

/**
 * The public coefficient object
 * 
 * @type {Object}
 */
Polynomial.prototype = {
  'coeff': null,
  /**
   * GCD over a field using Euclidean algorithm; returns monic gcd
   * @param {String|Object} x
   * @returns {Polynomial}
   */
  'gcd': function (x) {

    let a = normalize(clone(this['coeff']));
    let b = normalize(parse(x));

    if (isZeroPoly(a)) return new Polynomial(makeMonic(b, lcExp(b)));
    if (isZeroPoly(b)) return new Polynomial(makeMonic(a, lcExp(a)));

    while (!isZeroPoly(b)) {

      const r = normalize(clone(a));

      div(r, b); // r becomes remainder

      a = b;
      b = r;
    }

    return new Polynomial(makeMonic(a, lcExp(a)));
  },

  /**
   * Negate all coefficients of the polynomial
   * 
   * @returns {Polynomial}
   */
  'neg': function () {

    const ret = {};
    const poly = this['coeff'];

    for (const i in poly) {
      ret[i] = FIELD['mul'](poly[i], -1);
    }
    return new Polynomial(ret);
  },

  /**
   * Return the 'reciprocal polynomial', where the coefficients
   * appear in opposite order; i.e. a[i] -> a[n-i].
   * See e.g. https://en.wikipedia.org/wiki/Reciprocal_polynomial
   *
   * @returns {Polynomial}
   */
  'reciprocal': function () {

    const ret = {};
    const poly = this['coeff'];
    const n = degree(poly);

    for (const i in poly) {
      ret[n - i] = poly[i];
    }
    return new Polynomial(ret);
  },

  /**
   * Numerically evaluate the polynomial at a specific point x by
   * using Horner's method.
   * See e.g. https://en.wikipedia.org/wiki/Horner%27s_method
   *
   * @param {number} x The point where to evaluate this polynomial
   * @returns {number} The value P(x)
   */
  'eval': function (x) {

    const poly = this['coeff'];
    const n = degree(poly);

    if (n < 0) {
      return 0;
    }

    let ret = poly[n];

    for (let i = n - 1; i >= 0; i--) {
      ret = FIELD['mul'](ret, x);
      if (!FIELD['empty'](poly[i])) {
        ret = FIELD['add'](ret, poly[i]);
      }
    }
    return ret;
  },

  /**
   * Gets the leading coefficient
   * 
   * @returns {Polynomial}
   */
  'lc': function () {

    const e = lcExp(this['coeff']);
    return e === null ? undefined : this['coeff'][e];
  },

  /**
   * Gets the leading monomial
   * 
   * @returns {Polynomial}
   */
  'lm': function () {

    const e = lcExp(this['coeff']);
    const res = {};
    if (e !== null)
      res[e] = this['coeff'][e];
    return new Polynomial(res);
  },

  /**
   * Divide all coefficients of f by lc(f)
   * 
   * @returns {Polynomial}
   */
  'monic': function () {

    return new Polynomial(makeMonic(clone(this['coeff']), lcExp(this['coeff'])));
  },

  /**
   * Calculates the sum of two polynomials
   * 
   * @param {String|Object} x The summand polynomial
   * @returns {Polynomial}
   */
  'add': function (x) {

    const para = parse(x);

    const ret = {};
    const poly = this['coeff'];

    const keys = keyUnion(para, poly);

    for (const i in keys) {
      ret[i] = FIELD['add'](poly[i] || 0, para[i] || 0);
    }
    return new Polynomial(ret);
  },

  /**
   * Calculates the difference of two polynomials
   * 
   * @param {String|Object} x The subtrahend polynomial
   * @returns {Polynomial}
   */
  'sub': function (x) {

    const para = parse(x);

    const ret = {};
    const poly = this['coeff'];

    const keys = keyUnion(para, poly);

    for (const i in keys) {
      ret[i] = FIELD['sub'](poly[i] || 0, para[i] || 0);
    }
    return new Polynomial(ret);
  },

  /**
   * Calculates the product of two polynomials
   * 
   * @param {String|Object} x The minuend polynomial
   * @returns {Polynomial}
   */
  'mul': function (x) {

    const para = parse(x);

    const ret = {};
    const poly = this['coeff'];

    for (let i in para) {

      i = +i;

      for (let j in poly) {

        j = +j;

        ret[i + j] = FIELD['add'](ret[i + j] || 0, FIELD['mul'](para[i] || 0, poly[j] || 0));
      }
    }
    return new Polynomial(normalize(ret));
  },

  /**
   * Calculates the product of the two parameters and adds it to the current number (linear combination)
   * 
   * @param {String|Object} x The first factor polynomial
   * @param {String|Object} y The second factor polynomial
   * @returns {Polynomial}
   */
  'addmul': function (x, y) {

    const _x = parse(x);
    const _y = parse(y);

    const res = {};
    for (let i in _x) {

      i = +i;

      for (let j in _y) {
        j = +j;

        res[i + j] = FIELD['add'](res[i + j] || 0, FIELD['mul'](_x[i] || 0, _y[j] || 0));
      }
    }
    return this['add'](normalize(res));
  },

  /**
   * Calculates the quotient of two polynomials
   * 
   * @param {String|Object} x The denominator polynomial
   * @returns {Polynomial}
   */
  'div': function (x) {

    return new Polynomial(div(clone(this['coeff']), parse(x)));
  },

  /**
   * Calculates the pow of a polynomial to the exponent natural number e
   * 
   * @returns {Polynomial}
   */
  'pow': function (e) {

    if (isNaN(e) || e < 0 || e % 1) { // Only integer exponents
      throw new Error("Invalid");
    }

    let res = new Polynomial(1);
    let tmp = this;

    while (e > 0) {

      if (e & 1) {
        res = res['mul'](tmp);
      }
      tmp = tmp['mul'](tmp);
      e >>= 1;
    }
    return res;
  },

  /**
   * Calculates the modulo of a polynomial to another
   * 
   * @param {String|Object} x The second poly
   * @returns {Polynomial}
   */
  'mod': function (x) {

    const mod = clone(this['coeff']);

    div(mod, parse(x));

    return new Polynomial(mod);
  },

  /**
   * Calculates the nth derivative of the polynomial
   * 
   * @param {number} n The nth derivative
   * @returns {Polynomial}
   */
  'derive': function (n) {

    if (n === undefined) {
      n = 1;
    } else if (n < 0) {
      return null;
    }

    const poly = this['coeff'];
    const ret = {};

    for (const i in poly) {

      if (+i >= n)
        ret[i - n] = FIELD['mul'](poly[i] || 0, fallingFactorial(+i, n));
    }
    return new Polynomial(ret);
  },

  /**
   * Calculates the nth integral of the polynomial
   * 
   * @param {number} n The nth integral
   * @returns {Polynomial}
   */
  'integrate': function (n) {

    if (n === undefined) {
      n = 1;
    } else if (n < 0) {
      return null;
    }

    const poly = this['coeff'];
    const ret = {};

    for (const i in poly) {
      // divide by (i+1)(i+2)...(i+n)
      ret[+i + n] = FIELD['div'](poly[i] || 0, risingFactorial(+i, n));
    }
    return new Polynomial(ret);
  },

  /**
   * Formats the polynomial as a string
   * 
   * @returns {string} The polynomial string
   */
  'toString': toString("toString"),

  /**
   * Formats the polynomial as a latex representation
   * 
   * @returns {string} The polynomial latex string
   */
  'toLatex': toString("toLatex"),

  /**
   * Returns the actual polynomial in horner scheme
   * 
   * @returns {string}
   */
  'toHorner': function () {

    const poly = this['coeff'];
    const keys = [];
    for (const i in poly) {
      if (!FIELD.empty(poly[i]))
        keys.push(+i);
    }

    if (keys.length === 0)
      return "0";

    keys.sort(function (a, b) {
      return a - b;
    });

    // TODO: DRY, Combine with toString function
    function valToString(val, hasSign) {

      let str = "";

      if (Complex && val instanceof Complex) {

        if (val['im'] === 0) {

          if (val['re'] > 0 && hasSign) {
            str += "+";
          }
          str += val['re'];

        } else if (val['re'] === 0) {

          if (val['im'] === -1) {
            str += "-";
          } else if (val['im'] !== 1) {

            if (val['im'] > 0 && hasSign) {
              str += "+";
            }
            str += val['im'];
          } else {
            if (val['im'] > 0 && hasSign) {
              str += "+";
            }
          }
          str += "i";

        } else {

          if (hasSign) {
            str += "+";
          }

          str += "(";
          str += val.toString();
          str += ")";
        }

        return str;

      } else {

        if (val > 0 && hasSign) {
          str += "+";
        }
        str += val.toString();
      }
      return str;
    }

    function rec(keys, pos) {

      const ndx = keys.length - pos - 1;
      const exp = keys[ndx] - (keys[ndx - 1] || 0);
      let str1 = "";
      let str2 = "";

      if (exp > 0) {
        str1 = "x";

        if (exp > 1) {
          str1 += "^" + exp;
        }
      }

      if (ndx > 0)
        str1 += valToString(poly[keys[ndx - 1]], true);

      if (pos === 0) {
        return valToString(poly[keys[ndx]], false) + str1;
      }

      if (ndx >= 0 && keys[ndx])
        str2 += "(";

      str2 += rec(keys, pos - 1);

      if (ndx >= 0 && keys[ndx])
        str2 += ")";

      str2 += str1;

      return str2;
    }
    return rec(keys, keys.length - 1);
  },

  /**
   * Clones the actual object
   * 
   * @returns {Polynomial}
   */
  'clone': function () {
    return new Polynomial(this);
  },

  /**
   * Returns the degree of the polynomial
   * 
   * @returns {number}
   */
  'degree': function () {

    return degree(this['coeff']);
  }
};

/**
 * Form a (monic) polynomial out of an array of roots
 *
 * @param {Array<number>} roots - Array of roots
 * @returns {Polynomial} The monic polynomial with those roots
 */
Polynomial['fromRoots'] = function (roots) {

  const n = roots.length;

  const zero = FIELD['parse'](0);

  const nonZeroRoots = roots.filter(root => !(FIELD['equals'](root, zero)));
  const numZeros = n - nonZeroRoots.length;

  // First we construct the depressed polynomial with a recursive
  // strategy (this minimizes the number of multiplications)
  const pOne = new Polynomial(FIELD['parse'](1));

  function productHelper(r) {
    switch (r.length) {
      case 0:
        return pOne;
      case 1:
        return new Polynomial([FIELD['mul'](r[0], -1), 1]);
      default: // recurse
        const nLeft = Math.floor(r.length / 2);
        const left = r.slice(0, nLeft);
        const right = r.slice(nLeft, r.length);
        return productHelper(left)['mul'](productHelper(right));
    }
  }

  const dep = productHelper(nonZeroRoots);

  // Now raise the order by including numZeros zeros
  const dcoeff = dep['coeff'];
  const coeff = {};

  for (let i in dcoeff) {
    coeff[numZeros + parseInt(i, 10)] = dcoeff[i];
  }

  return new Polynomial(coeff);
};

/** 
 * True if polynomial map has no non-empty terms
 */
function isZeroPoly(r) {

  return degree(r) < 0;
}

/**
 * Set the field globally
 * 
 * @param {string|Object} field One of: C (complex), H (quaternion), Q (rational), R (real), Z<n>, or a custom field object
 */
Polynomial['setField'] = function (field) {

  // Fields with the same common API
  const F = {
    "Q": Fraction,
    "C": Complex,
    "H": Quaternion
  }[field];

  if (F !== undefined) {

    FIELD = {
      'add': function (a, b) { return new F(a)['add'](b); },
      'sub': function (a, b) { return new F(a)['sub'](b); },
      'mul': function (a, b) { return new F(a)['mul'](b); },
      'div': function (a, b) { return new F(a)['div'](b); },
      'parse': function (x) { return new F(x); },
      'empty': function (x) { return new F(x)['equals'](0); },
      'pow': function (a, b) { return new F(a)['pow'](b); },
      'equals': function (a, b) { return new F(a)['equals'](b); }
    };

  } else if (!field || field === 'R') {

    FIELD = FIELDS['R'];

  } else if (typeof field === 'object') {

    FIELD = field;

  } else if (field.charAt(0) === 'Z') {

    const N = +field.slice(1);

    FIELD = {// Test in Z_n
      'add': function (a, b) { return mod(a + b, N); },
      'sub': function (a, b) { return mod(a - b, N); },
      'mul': function (a, b) { return mod(a * b, N); },
      'div': function (a, b) { return mod(a * modinv(b, N), N); },
      'parse': function (x) { return parseInt(x, 10); },
      'empty': function (x) { return x === undefined || x === 0; },
      "pow": function (a, b) {

        let r = 1;
        while (b > 0) {

          if (b & 1) {
            r = mod(r * a, N);
          }
          a = mod(a * a, N);
          b >>= 1;
        }
        return r;
      },
      'equals': function (a, b) { return a == b; }
    };
  }
};
