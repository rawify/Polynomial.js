
declare module 'Polynomial';

class Polynomial {
  
  coeff: Record<string, any>;

  constructor(x: string | object);

  /**
   * Calculates the gcd of two polynomials
   * @param x The denominator polynomial
   * @returns Polynomial
   */
  gcd(x: string | object): Polynomial;

  /**
   * Negate all coefficients of the polynomial
   * @returns Polynomial
   */
  neg(): Polynomial;

  /**
   * Return the reciprocal polynomial, where the coefficients appear in opposite order
   * @returns Polynomial
   */
  reciprocal(): Polynomial;

  /**
   * Numerically evaluate the polynomial at a specific point x using Horner's method
   * @param x The point where to evaluate this polynomial
   * @returns number The value P(x)
   */
  eval(x: number): number;

  /**
   * Gets the leading coefficient
   * @returns Polynomial
   */
  lc(): Polynomial;

  /**
   * Gets the leading monomial
   * @returns Polynomial
   */
  lm(): Polynomial;

  /**
   * Divide all coefficients by lc(f)
   * @returns Polynomial
   */
  monic(): Polynomial;

  /**
   * Calculates the sum of two polynomials
   * @param x The summand polynomial
   * @returns Polynomial
   */
  add(x: string | object): Polynomial;

  /**
   * Calculates the difference of two polynomials
   * @param x The subtrahend polynomial
   * @returns Polynomial
   */
  sub(x: string | object): Polynomial;

  /**
   * Calculates the product of two polynomials
   * @param x The minuend polynomial
   * @returns Polynomial
   */
  mul(x: string | object): Polynomial;

  /**
   * Calculates the product of the two parameters and adds it to the current number (linear combination)
   * @param x The first factor polynomial
   * @param y The second factor polynomial
   * @returns Polynomial
   */
  addmul(x: string | object, y: string | object): Polynomial;

  /**
   * Calculates the quotient of two polynomials
   * @param x The denominator polynomial
   * @returns Polynomial
   */
  div(x: string | object): Polynomial;

  /**
   * Calculates the power of a polynomial to the exponent e
   * @param e The exponent
   * @returns Polynomial
   */
  pow(e: number): Polynomial;

  /**
   * Calculates the modulo of a polynomial to another
   * @param x The second polynomial
   * @returns Polynomial
   */
  mod(x: string | object): Polynomial;

  /**
   * Calculates the nth derivative of the polynomial
   * @param n The nth derivative
   * @returns Polynomial
   */
  derive(n: number): Polynomial;

  /**
   * Calculates the nth integral of the polynomial
   * @param n The nth integral
   * @returns Polynomial
   */
  integrate(n: number): Polynomial;

  /**
   * Formats the polynomial as a string
   * @returns string The polynomial string
   */
  toString(): string;

  /**
   * Formats the polynomial as a LaTeX representation
   * @returns string The polynomial LaTeX string
   */
  toLatex(): string;

  /**
   * Returns the actual polynomial in Horner scheme
   * @returns string
   */
  toHorner(): string;

  /**
   * Clones the polynomial object
   * @returns Polynomial
   */
  clone(): Polynomial;

  /**
   * Returns the degree of the polynomial
   * @returns number
   */
  degree(): number;

  /**
   * Set the field globally
   * @param field One of: C (complex), H (quaternion), Q (rational), R (real), or an object with methods for the field
   */
  static setField(field: string | object): void;

  /**
   * Form a (monic) polynomial out of an array of roots
   * @param roots Array of roots
   * @returns Polynomial The monic polynomial with those roots
   */
  static fromRoots(roots: number[]): Polynomial;
}

export { Polynomial as default, Polynomial };
