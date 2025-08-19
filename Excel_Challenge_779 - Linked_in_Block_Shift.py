# --- 1) CORE FUNCTIONS ---
# --- Paste this into a =PY() cell to define the helpers (no output expected from this cell): ---
# --- Rotation helpers ---
def rotate_right(lst, k):
    n = len(lst)
    if n == 0: 
        return lst[:]
    k %= n
    return lst[-k:] + lst[:-k] if k else lst[:]

def rotate_down_inplace(mat, col, k):
    r = len(mat)
    if r == 0: 
        return
    k %= r
    if k == 0: 
        return
    col_vals = [mat[i][col] for i in range(r)]
    new_vals = col_vals[-k:] + col_vals[:-k]
    for i in range(r):
        mat[i][col] = new_vals[i]

# --- One Set (rows right, then columns down) ---
def one_set_by_values(mat, row_shifts, col_shifts):
    R = len(mat)
    C = len(mat[0]) if R else 0
    out = [row[:] for row in mat]

    # Row right-rotations
    for r in range(R):
        k = row_shifts[r] % C if C else 0
        out[r] = rotate_right(out[r], k)

    # Column down-rotations
    for c in range(C):
        k = col_shifts[c] % R if R else 0
        rotate_down_inplace(out, c, k)

    return out

# --- Build permutation for one Set on row-major indices 0..R*C-1 ---
def perm_from_set(R, C, row_shifts, col_shifts):
    mat_idx = [[r*C + c for c in range(C)] for r in range(R)]
    moved = one_set_by_values(mat_idx, row_shifts, col_shifts)
    p = [0]*(R*C)   # p[i] = new index of item that was at i
    for r in range(R):
        for c in range(C):
            old_idx = moved[r][c]
            new_idx = r*C + c
            p[old_idx] = new_idx
    return p

# --- Permutation utilities ---
def perm_compose(p, q):
    # composition: first p, then q  -> result[i] = q[p[i]]
    return [q[p[i]] for i in range(len(p))]

def perm_power(p, k):
    # fast exponentiation of permutation p^k
    n = len(p)
    # identity
    result = list(range(n))
    base = p[:]
    while k > 0:
        if k & 1:
            result = perm_compose(result, base)
        base = perm_compose(base, base)
        k >>= 1
    return result

def cycles_and_order(p):
    n = len(p)
    seen = [False]*n
    lens = []
    for i in range(n):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L > 1:
                lens.append(L)
    order = math.lcm(*lens) if lens else 1
    return lens, order

# --- Apply a permutation to a matrix (row-major) ---
def apply_perm_to_matrix(mat, p):
    R = len(mat)
    C = len(mat[0]) if R else 0
    flat = sum(mat, [])
    new_flat = [None]*(R*C)
    for i, j in enumerate(p):         # item at i moves to j
        new_flat[j] = flat[i]
    # reshape
    out = [new_flat[r*C:(r+1)*C] for r in range(R)]
    return out
# --- Convenience: compute minimal equivalent sets and resulting grid ---
def compute_min_sets_and_grid(grid, row_shifts, col_shifts, target_sets):
    R = len(grid)
    C = len(grid[0]) if R else 0

    p = perm_from_set(R, C, row_shifts, col_shifts)
    _, order = cycles_and_order(p)

    kmin = target_sets % order
    pk = perm_power(p, kmin)
    out_grid = apply_perm_to_matrix(grid, pk)

    return kmin, out_grid, order
'Perform data manipulation'


# --- 2) Read from Excel, compute, and return results ---
# --- Paste this into another =PY() cell. It reads your sheet ranges (matching your layout), computes the minimum sets equivalent to 50,000, and returns the resulting grid as a DataFrame with headers. You can change target_sets to any number you need. ---

import math
# --- Read inputs from the current sheet ---
grid = np.array(xl("B2:F4")).tolist()          # letters grid (3x5)
row_shifts = [int(x) for x in np.array(xl("A2:A4")).flatten().tolist()]   # e.g., [1,4,3]
col_shifts = [int(x) for x in np.array(xl("B1:F1")).flatten().tolist()]   # e.g., [2,3,1,5,4]

# --- How many Sets do you want to match? ---
target_sets = 50000   # change as needed

# --- Compute ---
kmin, out_grid, order = compute_min_sets_and_grid(grid, row_shifts, col_shifts, target_sets)

# Show kmin and order back into Excel in convenient cells (optional prints for your reference)
print(f"Permutation order: {order}")
print(f"Minimum equivalent Sets for {target_sets}: {kmin}")

# Return the resulting grid
out_grid
