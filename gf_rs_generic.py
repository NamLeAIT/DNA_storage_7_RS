from __future__ import annotations

class GF2m:
    def __init__(self, m: int, primitive_poly: int):
        self.m = int(m)
        self.primitive_poly = int(primitive_poly)
        self.size = 1 << self.m
        self.mask = self.size - 1
        self.order = self.size - 1
        self.exp = [0] * (self.order * 2)
        self.log = [0] * self.size

        x = 1
        for i in range(self.order):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & self.size:
                x ^= self.primitive_poly
            x &= self.mask

        for i in range(self.order, self.order * 2):
            self.exp[i] = self.exp[i - self.order]

    def add(self, x: int, y: int) -> int:
        return (int(x) ^ int(y)) & self.mask

    sub = add

    def mul(self, x: int, y: int) -> int:
        x &= self.mask
        y &= self.mask
        if x == 0 or y == 0:
            return 0
        return self.exp[self.log[x] + self.log[y]]

    def div(self, x: int, y: int) -> int:
        x &= self.mask
        y &= self.mask
        if y == 0:
            raise ZeroDivisionError("GF division by zero")
        if x == 0:
            return 0
        return self.exp[(self.log[x] - self.log[y]) % self.order]

    def pow(self, x: int, power: int) -> int:
        x &= self.mask
        if x == 0:
            return 0
        return self.exp[(self.log[x] * int(power)) % self.order]

    def inverse(self, x: int) -> int:
        x &= self.mask
        if x == 0:
            raise ZeroDivisionError("GF inverse of zero")
        return self.exp[(self.order - self.log[x]) % self.order]

    def poly_scale(self, p: list[int], x: int) -> list[int]:
        return [self.mul(coef, x) for coef in p]

    def poly_add(self, p: list[int], q: list[int]) -> list[int]:
        r = [0] * max(len(p), len(q))
        for i in range(len(p)):
            r[i + len(r) - len(p)] = p[i] & self.mask
        for i in range(len(q)):
            r[i + len(r) - len(q)] ^= q[i] & self.mask
        return r

    def poly_mul(self, p: list[int], q: list[int]) -> list[int]:
        r = [0] * (len(p) + len(q) - 1)
        for j, qj in enumerate(q):
            for i, pi in enumerate(p):
                r[i + j] ^= self.mul(pi, qj)
        return r

    def poly_div(self, dividend: list[int], divisor: list[int]) -> tuple[list[int], list[int]]:
        msg_out = list(dividend)
        for i in range(len(dividend) - len(divisor) + 1):
            coef = msg_out[i]
            if coef:
                for j in range(1, len(divisor)):
                    if divisor[j]:
                        msg_out[i + j] ^= self.mul(divisor[j], coef)
        sep = -(len(divisor) - 1)
        return msg_out[:sep], msg_out[sep:]

    def poly_eval(self, poly: list[int], x: int) -> int:
        y = poly[0]
        for i in range(1, len(poly)):
            y = self.mul(y, x) ^ poly[i]
        return y

    def generator_poly(self, nsym: int) -> list[int]:
        g = [1]
        for i in range(nsym):
            g = self.poly_mul(g, [1, self.pow(2, i)])
        return g

    def rs_encode_msg(self, msg: list[int], nsym: int) -> list[int]:
        msg = [int(x) & self.mask for x in msg]
        if len(msg) + nsym > self.order:
            raise ValueError("RS message too long for this field.")
        gen = self.generator_poly(nsym)
        msg_out = msg + [0] * nsym
        for i in range(len(msg)):
            coef = msg_out[i]
            if coef:
                for j in range(1, len(gen)):
                    msg_out[i + j] ^= self.mul(gen[j], coef)
        return msg + msg_out[-nsym:]

    def rs_calc_syndromes(self, msg: list[int], nsym: int) -> list[int]:
        return [0] + [self.poly_eval(msg, self.pow(2, i)) for i in range(nsym)]

    def rs_find_error_locator(self, synd: list[int], nsym: int, erase_count: int = 0) -> list[int]:
        err_loc = [1]
        old_loc = [1]
        synd_shift = len(synd) - nsym
        for i in range(nsym - erase_count):
            k = i + synd_shift
            delta = synd[k]
            for j in range(1, len(err_loc)):
                delta ^= self.mul(err_loc[-(j + 1)], synd[k - j])
            old_loc.append(0)
            if delta:
                if len(old_loc) > len(err_loc):
                    new_loc = self.poly_scale(old_loc, delta)
                    old_loc = self.poly_scale(err_loc, self.inverse(delta))
                    err_loc = new_loc
                err_loc = self.poly_add(err_loc, self.poly_scale(old_loc, delta))

        while len(err_loc) and err_loc[0] == 0:
            del err_loc[0]

        errs = len(err_loc) - 1
        if errs * 2 > nsym:
            raise ValueError("Too many errors to correct.")
        return err_loc

    def rs_find_errors(self, err_loc: list[int], nmess: int) -> list[int]:
        errs = len(err_loc) - 1
        err_pos = []
        for i in range(nmess):
            if self.poly_eval(err_loc, self.pow(2, i)) == 0:
                err_pos.append(nmess - 1 - i)
        if len(err_pos) != errs:
            raise ValueError("Could not locate errors.")
        return err_pos

    def rs_find_errata_locator(self, e_pos: list[int]) -> list[int]:
        e_loc = [1]
        for i in e_pos:
            e_loc = self.poly_mul(e_loc, self.poly_add([1], [self.pow(2, i), 0]))
        return e_loc

    def rs_find_error_evaluator(self, synd: list[int], err_loc: list[int], nsym: int) -> list[int]:
        _, remainder = self.poly_div(self.poly_mul(synd, err_loc), [1] + [0] * (nsym + 1))
        return remainder

    def rs_correct_errata(self, msg: list[int], synd: list[int], err_pos: list[int]) -> list[int]:
        coef_pos = [len(msg) - 1 - p for p in err_pos]
        err_loc = self.rs_find_errata_locator(coef_pos)
        err_eval = self.rs_find_error_evaluator(synd[::-1], err_loc, len(err_loc) - 1)[::-1]
        xs = [self.pow(2, -(self.order - pos)) for pos in coef_pos]
        e = [0] * len(msg)

        for i, xi in enumerate(xs):
            xi_inv = self.inverse(xi)
            err_loc_prime = 1
            for j, xj in enumerate(xs):
                if j != i:
                    err_loc_prime = self.mul(err_loc_prime, 1 ^ self.mul(xi_inv, xj))
            y = self.poly_eval(err_eval[::-1], xi_inv)
            y = self.mul(xi, y)
            e[err_pos[i]] = self.div(y, err_loc_prime)

        return self.poly_add(msg, e)

    def rs_correct_msg(self, msg_in: list[int], nsym: int) -> list[int]:
        msg = [int(x) & self.mask for x in msg_in]
        if len(msg) > self.order:
            raise ValueError("RS message too long.")

        synd = self.rs_calc_syndromes(msg, nsym)
        if max(synd) == 0:
            return msg[:-nsym]

        err_loc = self.rs_find_error_locator(synd, nsym)
        err_pos = self.rs_find_errors(err_loc[::-1], len(msg))
        corrected = self.rs_correct_errata(msg, synd, err_pos)

        if max(self.rs_calc_syndromes(corrected, nsym)) != 0:
            raise ValueError("Reed-Solomon correction failed.")

        return corrected[:-nsym]

    def lagrange_eval(self, xs: list[int], ys: list[int], x: int) -> int:
        total = 0
        for i, xi in enumerate(xs):
            num = 1
            den = 1
            for j, xj in enumerate(xs):
                if i == j:
                    continue
                num = self.mul(num, self.sub(x, xj))
                den = self.mul(den, self.sub(xi, xj))
            total ^= self.mul(ys[i], self.div(num, den))
        return total
