import type ZqField from "./f1field";
import * as Scalar from "./scalar";
// Check here: https://eprint.iacr.org/2012/685.pdf

export default function buildSqrt(F: ZqField): void {
  if (F.m % 2 == 1) {
    if (Scalar.eq(Scalar.mod(F.p, 4), 1)) {
      if (Scalar.eq(Scalar.mod(F.p, 8), 1)) {
        if (Scalar.eq(Scalar.mod(F.p, 16), 1)) {
          // alg7_muller(F);
          alg5_tonelliShanks(F);
        } else if (Scalar.eq(Scalar.mod(F.p, 16), 9)) {
          alg4_kong(F);
        } else {
          throw new Error("Field withot sqrt");
        }
      } else if (Scalar.eq(Scalar.mod(F.p, 8), 5)) {
        alg3_atkin(F);
      } else {
        throw new Error("Field withot sqrt");
      }
    } else if (Scalar.eq(Scalar.mod(F.p, 4), 3)) {
      alg2_shanks(F);
    }
  } else {
    const pm2mod4 = Scalar.mod(Scalar.pow(F.p, F.m / 2), 4);
    if (pm2mod4 == 1n) {
      alg10_adj(F);
    } else if (pm2mod4 == 3n) {
      alg9_adj(F);
    } else {
      alg8_complex(F);
    }
  }
}

function alg5_tonelliShanks(F: ZqField): void {
  F.sqrt_q = Scalar.pow(F.p, F.m);

  F.sqrt_s = 0;
  F.sqrt_t = Scalar.sub(F.sqrt_q, 1);

  while (!Scalar.isOdd(F.sqrt_t)) {
    F.sqrt_s = F.sqrt_s + 1;
    F.sqrt_t = Scalar.div(F.sqrt_t, 2);
  }

  let c0 = F.one;

  while (F.eq(c0, F.one)) {
    const c = F.random();
    F.sqrt_z = F.pow(c, F.sqrt_t);
    c0 = F.pow(F.sqrt_z, 2 ** (F.sqrt_s - 1));
  }

  F.sqrt_tm1d2 = Scalar.div(Scalar.sub(F.sqrt_t, 1), 2);

  F.sqrt = function (this: ZqField, a: bigint): bigint | null {
    if (this.isZero(a)) return this.zero;
    let w = this.pow(a, this.sqrt_tm1d2!);
    const a0 = this.pow(this.mul(this.square(w), a), 2 ** (this.sqrt_s! - 1));
    if (this.eq(a0, this.negone)) return null;

    let v = this.sqrt_s!;
    let x = this.mul(a, w);
    let b = this.mul(x, w);
    let z = this.sqrt_z!;
    while (!this.eq(b, this.one)) {
      let b2k = this.square(b);
      let k = 1;
      while (!this.eq(b2k, this.one)) {
        b2k = this.square(b2k);
        k++;
      }

      w = z;
      for (let i = 0; i < v - k - 1; i++) {
        w = this.square(w);
      }
      z = this.square(w);
      b = this.mul(b, z);
      x = this.mul(x, w);
      v = k;
    }
    return this.geq(x, this.zero) ? x : this.neg(x);
  };
}

function alg4_kong(F: ZqField): void {
  F.sqrt = function (): bigint | null {
    throw new Error("Sqrt alg 4 not implemented");
  };
}

function alg3_atkin(F: ZqField): void {
  F.sqrt = function (): bigint | null {
    throw new Error("Sqrt alg 3 not implemented");
  };
}

function alg2_shanks(F: ZqField): void {
  F.sqrt_q = Scalar.pow(F.p, F.m);
  F.sqrt_e1 = Scalar.div(Scalar.sub(F.sqrt_q, 3), 4);

  F.sqrt = function (this: ZqField, a: bigint): bigint | null {
    if (this.isZero(a)) return this.zero;

    // Test that have solution
    const a1 = this.pow(a, this.sqrt_e1!);

    const a0 = this.mul(this.square(a1), a);

    if (this.eq(a0, this.negone)) return null;

    const x = this.mul(a1, a);

    return this.geq(x, this.zero) ? x : this.neg(x);
  };
}

function alg10_adj(F: ZqField): void {
  F.sqrt = function (): bigint | null {
    throw new Error("Sqrt alg 10 not implemented");
  };
}

function alg9_adj(F: ZqField): void {
  F.sqrt_q = Scalar.pow(F.p, F.m / 2);
  F.sqrt_e34 = Scalar.div(Scalar.sub(F.sqrt_q, 3), 4);
  F.sqrt_e12 = Scalar.div(Scalar.sub(F.sqrt_q, 1), 2);

  F.frobenius = function (n: number, x: bigint): bigint {
    if (n % 2 == 1) {
      return (F as any).conjugate(x);
    } else {
      return x;
    }
  };

  F.sqrt = function (this: ZqField, a: bigint): bigint | null {
    const a1 = this.pow(a, this.sqrt_e34!);
    const alfa = this.mul(this.square(a1), a);
    const a0 = this.mul(this.frobenius!(1, alfa), alfa);
    if (this.eq(a0, this.negone)) return null;
    const x0 = this.mul(a1, a);
    let x: bigint;
    if (this.eq(alfa, this.negone)) {
      x = this.mul(x0, [(this as any).F.zero, (this as any).F.one]);
    } else {
      const b = this.pow(this.add(this.one, alfa), this.sqrt_e12!);
      x = this.mul(b, x0);
    }
    return this.geq(x, this.zero) ? x : this.neg(x);
  };
}

function alg8_complex(F: ZqField): void {
  F.sqrt = function (): bigint | null {
    throw new Error("Sqrt alg 8 not implemented");
  };
}
