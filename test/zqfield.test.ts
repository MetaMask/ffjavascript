import { describe, expect, it } from "vitest";

import ZqField from "../src/f1field";
import * as Scalar from "../src/scalar";

const q = Scalar.fromString("21888242871839275222246405745257275088696311157297823662689037894645226208583")!;
const r = Scalar.fromString("21888242871839275222246405745257275088548364400416034343698204186575808495617")!;

describe("F1 testing", () => {
  it("Should compute euclidean", () => {
    const F = new ZqField(7);
    const res = F.inv(F.e(4));

    expect(F.eq(res, F.e(2))).toBe(true);
  });

  it("Should multiply and divide in F1", () => {
    const F = new ZqField(q);
    const a = F.e("1");
    const b = F.normalize(-3);
    const c = F.mul(a, b);
    const d = F.div(c, b);

    expect(F.eq(a, d)).toBe(true);
  });

  it("Should compute sqrts", () => {
    const F = new ZqField(q);
    const a = F.e("4");
    const b = F.sqrt!(a);

    expect(F.eq(F.e(0), F.sqrt!(F.e("0"))!)).toBe(true);
    expect(F.eq(b!, F.e("2"))).toBe(true);
    expect(F.sqrt!(F.nqr)).toBe(null);
  });

  it("Should compute sqrt of 100 random numbers", () => {
    const F = new ZqField(r);
    for (let j = 0; j < 100; j++) {
      const a = F.random();
      const s = F.sqrt!(a);
      if (s != null) {
        expect(F.eq(F.square(s), a)).toBe(true);
      }
    }
  });
});
