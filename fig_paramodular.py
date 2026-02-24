"""
fig_paramodular.py — Figures for Paper #15
"""
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for d in range(3, int(n**0.5)+1, 2):
        if n % d == 0: return False
    return True

def count_points(p, q, r):
    count = 1  # infinity
    for x in range(r):
        val = (x * (x*x - p*p) * (x*x - q*q)) % r
        if val == 0:
            count += 1
        elif pow(val, (r-1)//2, r) == 1:
            count += 2
    return count

# Collect traces for (p,q) = (3,7)
p0, q0 = 3, 7
bad = {2, 3, 5, 7}

data_r1 = []  # (r, a_r) for r ≡ 1 mod 4
data_r3 = []  # (r, a_r) for r ≡ 3 mod 4

for r in range(11, 3000):
    if not is_prime(r): continue
    if r in bad: continue
    nr = count_points(p0, q0, r)
    a_r = r + 1 - nr
    if r % 4 == 1:
        data_r1.append((r, a_r))
    else:
        data_r3.append((r, a_r))

print(f"Collected: {len(data_r1)} primes ≡ 1 (mod 4), {len(data_r3)} primes ≡ 3 (mod 4)")

# Also collect for multiple curves
curves = [(3,7),(7,23),(11,19),(13,17),(3,17)]
multi_data = {}
for p, q in curves:
    N = (p+q)//2; diff = abs(p-q)
    bad_set = {2}
    for n in [p, q, N, diff]:
        d = 2; temp = abs(n)
        while d*d <= temp:
            if temp % d == 0:
                bad_set.add(d)
                while temp % d == 0: temp //= d
            d += 1
        if temp > 1: bad_set.add(temp)
    
    traces_1 = []; traces_3 = []
    for r in range(11, 2000):
        if not is_prime(r) or r in bad_set: continue
        nr = count_points(p, q, r)
        a_r = r + 1 - nr
        if r % 4 == 1: traces_1.append(a_r)
        else: traces_3.append(a_r)
    multi_data[(p,q)] = (traces_1, traces_3)

# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Trace vanishing law
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel (a): a_r vs r, colored by r mod 4
ax = axes[0]
rs1 = [d[0] for d in data_r1]
as1 = [d[1] for d in data_r1]
rs3 = [d[0] for d in data_r3]
as3 = [d[1] for d in data_r3]

ax.scatter(rs3, as3, s=6, c='#CC3333', alpha=0.5, label=r'$r \equiv 3\ (\mathrm{mod}\ 4)$: $a_r = 0$ always', zorder=5)
ax.scatter(rs1, as1, s=6, c='#2255BB', alpha=0.3, label=r'$r \equiv 1\ (\mathrm{mod}\ 4)$: $a_r$ varies', zorder=4)
ax.axhline(y=0, color='gray', lw=0.5, alpha=0.5)
ax.set_xlabel(r'Good prime $r$', fontsize=12)
ax.set_ylabel(r'Frobenius trace $a_r$', fontsize=12)
ax.set_title(r'(a) Trace vanishing law for $(p,q) = (3,7)$', fontsize=12)
ax.legend(fontsize=8, loc='upper left')
ax.set_xlim(0, 3000)
ax.grid(True, alpha=0.1)

# Panel (b): Histogram of normalised traces
ax = axes[1]
norm_r1 = [a / (2*math.sqrt(r)) for r, a in data_r1 if a != 0]
norm_all = [a / (2*math.sqrt(r)) for r, a in data_r1]

# USp(4) Sato-Tate would have smooth distribution
# Our distribution is concentrated at 0 with tails
bins = np.linspace(-2, 2, 41)
ax.hist(norm_all, bins=bins, color='#2255BB', alpha=0.6, edgecolor='navy',
       label=r'$r \equiv 1\ (\mathrm{mod}\ 4)$', density=True)

# Add the r ≡ 3 spike at 0
# These are all exactly 0, so add a bar at center
r3_frac = len(data_r3) / (len(data_r1) + len(data_r3))
ax.axvline(x=0, color='#CC3333', lw=2, alpha=0.7,
          label=r'$r \equiv 3\ (\mathrm{mod}\ 4)$: all at $0$')
ax.annotate(f'{len(data_r3)} primes\n(all $a_r = 0$)', 
           xy=(0, 0), xytext=(0.8, 3.5),
           fontsize=9, color='#CC3333',
           arrowprops=dict(arrowstyle='->', color='#CC3333'))

ax.set_xlabel(r'Normalised trace $a_r / 2\sqrt{r}$', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title(r'(b) Trace distribution for $(p,q) = (3,7)$', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.1)

plt.tight_layout()
plt.savefig('/home/claude/paper15/figures/fig_trace_vanishing.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/paper15/figures/fig_trace_vanishing.png', dpi=200, bbox_inches='tight')
plt.close()
print("Figure 1 done.")


# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Cross-validation across 5 curves + local L-factors
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel (a): Vanishing rate by curve
ax = axes[0]
labels = []
zeros_r1 = []
zeros_r3 = []
total_r1 = []
total_r3 = []

for (p, q), (t1, t3) in multi_data.items():
    labels.append(f'({p},{q})')
    z1 = sum(1 for t in t1 if t == 0)
    zeros_r1.append(z1 / len(t1) * 100 if t1 else 0)
    zeros_r3.append(sum(1 for t in t3 if t == 0) / len(t3) * 100 if t3 else 0)
    total_r1.append(len(t1))
    total_r3.append(len(t3))

x = np.arange(len(labels))
width = 0.35
bars1 = ax.bar(x - width/2, zeros_r1, width, color='#2255BB', alpha=0.7,
              label=r'$r \equiv 1\ (\mathrm{mod}\ 4)$')
bars3 = ax.bar(x + width/2, zeros_r3, width, color='#CC3333', alpha=0.7,
              label=r'$r \equiv 3\ (\mathrm{mod}\ 4)$')

ax.set_ylabel(r'$\%$ of primes with $a_r = 0$', fontsize=12)
ax.set_title(r'(a) Vanishing rate across 5 curves', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylim(0, 110)
ax.axhline(y=100, color='#CC3333', lw=1, ls=':', alpha=0.5)
ax.axhline(y=50, color='#2255BB', lw=1, ls=':', alpha=0.5)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.1, axis='y')

# Panel (b): Local L-factor type summary  
ax = axes[1]
ax.axis('off')

# Create a nice summary table
table_data = [
    ['Case', 'Condition', 'Reduction', 'Local $L$-factor'],
    ['I/II', '$r = p$ or $q$', '1 cusp ($A_2$)', '$(1 - a_E r^{-s} + r^{1-2s})^{-1}$'],
    ['III/IV', '$r \\mid N$ or $r \\mid (p{-}q)$', '2 nodes ($A_1$)', '$(1 - r^{-s})^{-2}$'],
    ['Good, $r \\equiv 3$', '$r \\nmid \\Delta$', 'Smooth', '$P_r(r^{-s})^{-1}$, $a_1 = 0$'],
    ['Good, $r \\equiv 1$', '$r \\nmid \\Delta$', 'Smooth', '$P_r(r^{-s})^{-1}$, generic'],
]

# Draw as text
y_start = 0.95
ax.text(0.5, y_start, 'Local $L$-factor classification', fontsize=13, 
       fontweight='bold', ha='center', va='top', transform=ax.transAxes)

colors_row = ['#f0f0f0', '#ffffff', '#f0f0f0', '#ffffff']
headers = ['Case', 'Condition', 'Reduction', 'Local factor']
col_x = [0.02, 0.22, 0.50, 0.72]

for j, h in enumerate(headers):
    ax.text(col_x[j], 0.82, h, fontsize=10, fontweight='bold',
           transform=ax.transAxes, va='top')

row_data = [
    ['I / II', 'r = p or q', '1 cusp (A₂)', '(1-aₑr⁻ˢ+r¹⁻²ˢ)⁻¹'],
    ['III / IV', 'r | N or r | (p-q)', '2 nodes (A₁)', '(1-r⁻ˢ)⁻²'],
    ['Good, r≡3(4)', 'r ∤ Δ', 'Smooth', 'a₁ = 0 (proved)'],
    ['Good, r≡1(4)', 'r ∤ Δ', 'Smooth', 'a₁ generic'],
]

for i, row in enumerate(row_data):
    y = 0.68 - i * 0.14
    bg_color = '#E8E8FF' if i < 2 else '#FFE8E8' if i == 2 else '#E8FFE8'
    ax.axhspan(y - 0.05, y + 0.08, xmin=0, xmax=1, 
              facecolor=bg_color, alpha=0.5, transform=ax.transAxes)
    for j, val in enumerate(row):
        ax.text(col_x[j], y, val, fontsize=9, transform=ax.transAxes, va='center')

ax.set_title('(b) Complete local arithmetic at all primes', fontsize=12)

plt.tight_layout()
plt.savefig('/home/claude/paper15/figures/fig_cross_validate.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/paper15/figures/fig_cross_validate.png', dpi=200, bbox_inches='tight')
plt.close()
print("Figure 2 done.")
