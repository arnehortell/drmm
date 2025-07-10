
import pandas as pd
from sympy import isprime as miller_rabin_test

def hybrid_drmm_with_miller_rabin(n):
    n_bin = bin(n)[2:]
    n_bits = list(map(int, reversed(n_bin)))
    result = None
    path_count = 0
    used_mr = False

    for a in range(3, int(n**0.5) + 1, 2):  # hoppa över jämna a
        path_count += 1

        if path_count == 5:
            if miller_rabin_test(n):
                used_mr = True
                return {
                    "binary": n_bin,
                    "is_prime": True,
                    "factors": None,
                    "paths_tested": path_count,
                    "used_miller_rabin": used_mr
                }

        if n % a != 0:
            continue

        b = n // a

        a_bin = list(map(int, reversed(bin(a)[2:])))
        b_bin = list(map(int, reversed(bin(b)[2:])))
        len_a = len(a_bin)
        len_b = len(b_bin)
        size = max(len_a, len_b)

        a_pad = a_bin + [0] * (size - len_a)
        b_pad = b_bin + [0] * (size - len_b)

        product_bits = [0] * (2 * size)
        for i in range(size):
            for j in range(size):
                product_bits[i + j] += a_pad[i] * b_pad[j]

        for i in range(len(product_bits)):
            if product_bits[i] >= 2:
                carry = product_bits[i] // 2
                product_bits[i] = product_bits[i] % 2
                if i + 1 < len(product_bits):
                    product_bits[i + 1] += carry
                else:
                    product_bits.append(carry)

        while len(product_bits) > 1 and product_bits[-1] == 0:
            product_bits.pop()

        if product_bits == n_bits:
            result = (a, b)
            break

    return {
        "binary": n_bin,
        "is_prime": result is None,
        "factors": result,
        "paths_tested": path_count,
        "used_miller_rabin": used_mr
    }

if __name__ == "__main__":
    results = []
    for n in range(1000001, 1001001, 2):  # bara udda tal
        res = hybrid_drmm_with_miller_rabin(n)
        res["decimal"] = n
        results.append(res)

    df = pd.DataFrame(results)
    df["factor_a"] = df["factors"].apply(lambda x: x[0] if x else None)
    df["factor_b"] = df["factors"].apply(lambda x: x[1] if x else None)

    num_mr_used = df["used_miller_rabin"].sum()
    num_total = len(df)
    num_primes = df["is_prime"].sum()
    avg_paths = df["paths_tested"].mean()
    avg_paths_mr = df[df["used_miller_rabin"]]["paths_tested"].mean()
    avg_paths_no_mr = df[~df["used_miller_rabin"]]["paths_tested"].mean()

    with open("hybrid_drmm_stats.txt", "w") as f:
        f.write(f"Totalt testade tal: {num_total}\n")
        f.write(f"Primtal identifierade: {num_primes}\n")
        f.write(f"Miller–Rabin användes: {num_mr_used}\n")
        f.write(f"Genomsnittliga paths (alla): {avg_paths:.2f}\n")
        f.write(f"Genomsnittliga paths (MR-fall): {avg_paths_mr:.2f}\n")
        f.write(f"Genomsnittliga paths (endast DRMM): {avg_paths_no_mr:.2f}\n")

    df.to_csv("hybrid_drmm_results.csv", index=False)
    print("Tabell sparad till hybrid_drmm_results.csv")
    print("Statistik sparad till hybrid_drmm_stats.txt")
