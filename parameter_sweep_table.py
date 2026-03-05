
# Summary table
import pandas as pd

df = pd.DataFrame(results)
df = df.sort_values('min_P2').reset_index(drop=True)

# ── Formatting helpers
def fmt_lambda(x):
    if x < 0.01:
        return f"{x:.4f}"
    elif x < 1:
        return f"{x:.3f}"
    elif x < 10:
        return f"{x:.2f}"
    else:
        return f"{x:.1f}"

def fmt_R(x):
    return f"{x:.1f}"

def fmt_P2(x):
    return f"{x:.6e}"

def fmt_rho(x):
    return f"{x:.4f}"

#Table (a): 30 configurations with smallest min P2
top30 = df.head(30).copy()
top30['lambda_fmt']     = top30['lambda'].map(fmt_lambda)
top30['R_fmt']          = top30['R'].map(fmt_R)
top30['min_P2_fmt']     = top30['min_P2'].map(fmt_P2)
top30['rho_at_min_fmt'] = top30['rho_at_min'].map(fmt_rho)
top30['P2 < 0']         = top30['P2_negative'].map(lambda b: '✓' if b else '✗')

display_cols = {
    'N'            : 'N',
    'lambda_fmt'   : 'λ',
    'R_fmt'        : 'R',
    'rho_at_min_fmt': 'ρ* at min',
    'min_P2_fmt'   : 'min P₂(ρ*)',
    'P2 < 0'       : 'P₂ < 0?',
}

print("=" * 70)
print("Table 1: 30 parameter configurations with the smallest min P₂(ρ*)")
print("         (sorted ascending; ✓ = P₂ strictly negative somewhere)")
print("=" * 70)
display(
    top30[list(display_cols.keys())]
    .rename(columns=display_cols)
    .style
    .set_caption("Minimum P₂(ρ*) over ρ* ∈ (0.05, 0.95) — 30 most extreme configurations")
    .set_properties(**{'text-align': 'right'})
    .apply(lambda col: [
        'background-color: #d4edda; color: #155724' if v == '✓' else ''
        for v in col
    ], subset=['P₂ < 0?'])
    .format(na_rep='—')
)

#Table (b): all configurations where P2 < 0
neg_df = df[df['P2_negative']].copy()

print()
print("=" * 70)
if neg_df.empty:
    print("Table 2: NO configurations found with P₂ < 0")
    print()
    print("Conclusion: P₂(ρ*) ≥ 0 for all (N, λ, R, ρ*) in the parameter space")
    print("explored. Informational dominance FAILS throughout — a more precise")
    print("prior is a strategic disadvantage at every symmetric equilibrium.")
else:
    print(f"Table 2: All {len(neg_df)} configurations where P₂ < 0 (informational dominance HOLDS)")
    print("=" * 70)
    neg_df['lambda_fmt']     = neg_df['lambda'].map(fmt_lambda)
    neg_df['R_fmt']          = neg_df['R'].map(fmt_R)
    neg_df['min_P2_fmt']     = neg_df['min_P2'].map(fmt_P2)
    neg_df['rho_at_min_fmt'] = neg_df['rho_at_min'].map(fmt_rho)

    display(
        neg_df[['N', 'lambda_fmt', 'R_fmt', 'rho_at_min_fmt', 'min_P2_fmt']]
        .rename(columns={
            'lambda_fmt'    : 'λ',
            'R_fmt'         : 'R',
            'rho_at_min_fmt': 'ρ* at min',
            'min_P2_fmt'    : 'min P₂(ρ*)',
        })
        .style
        .set_caption("Configurations where P₂(ρ*) < 0: informational dominance holds")
        .set_properties(**{'text-align': 'right'})
        .format(na_rep='—')
    )
print("=" * 70)

#Global extremes
print()
print("Global statistics across all valid configurations:")
print(f"  Overall minimum of min P₂ : {df['min_P2'].min():.6e}  "
      f"(N={df.loc[df['min_P2'].idxmin(),'N']}, "
      f"λ={fmt_lambda(df.loc[df['min_P2'].idxmin(),'lambda'])}, "
      f"R={fmt_R(df.loc[df['min_P2'].idxmin(),'R'])}, "
      f"ρ*={fmt_rho(df.loc[df['min_P2'].idxmin(),'rho_at_min'])})")
print(f"  Overall maximum of min P₂ : {df['min_P2'].max():.6e}")
print(f"  Fraction with P₂ < 0      : {df['P2_negative'].sum()} / {len(df)}  "
      f"({100*df['P2_negative'].mean():.2f}%)")
