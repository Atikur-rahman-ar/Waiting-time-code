#!/usr/bin/env python3
"""
waiting_time_analysis.py
========================
Single deposited numerical pipeline for the paper

  "The Waiting-Time Problem for Evolutionary Innovation:
   A Quantitative Synthesis of Population-Genetic Constraints on Macroevolution, 1980-2026"

Every waiting time, table cell, scenario estimate and appendix number quoted in the
manuscript is emitted by this script. Conventions are fixed once, here:

  * Ne        : effective population size in DIPLOID individuals (hominin central value 1e4)
  * 2*Ne      : number of gene copies producing mutations each generation
  * mu        : per-site mutation rate per generation (1.25e-8)
  * B         : time budget in generations (2.4e5 = 6 Myr at 25 yr/generation)
  * omega_k   : C_k / 4**k with C_k = number of functional k-tuples among the 4**k words
  * T_arr(k)  : 1 / (2*Ne * mu**k * omega_k)      [de novo, simultaneous, single specified locus]
  * P_seg(k)  : <= 2*Ne * mu**k * omega_k * 4*Ne  [segregating probability, neutral sojourn tau = 4Ne]

Published numbers quoted from other papers (Durrett & Schmidt 2008; Sanford et al. 2015;
Lynch 2005; Kimura & Ohta 1969) are reproduced here to the precision shown in those papers,
and are labelled QUOTED. All other numbers are COMPUTED by the formulas above.

Output: a plain-text log (numerical_outputs.txt) deposited with the manuscript as
Supplementary File S2.
"""
import math

# ----------------------------------------------------------------------
# Parameters (Table 1 of the manuscript)
# ----------------------------------------------------------------------
Ne      = 1.0e4          # diploid effective size (hominin, central)
Ncopies = 2.0 * Ne       # gene copies per generation  = 2e4
mu      = 1.25e-8        # per-site mutation rate per generation
B       = 2.4e5          # budget, generations (6 Myr / 25 yr)
B_gen   = 4.5e5          # generous budget (9 Myr / 20 yr)
YR      = 25.0           # years per generation
log10   = math.log10

L = []                   # output log


def p(line=""):
    L.append(str(line))
    print(line)


def sec(title):
    p()
    p("=" * 78)
    p(title)
    p("=" * 78)


def T_arr(k, C, N_eff=Ne, mu_=mu):
    """Eq. (3)-(4): first-arrival time, de novo simultaneous, single specified locus."""
    omega = C / (4.0 ** k)
    return 1.0 / (2.0 * N_eff * mu_ ** k * omega)


def lg(x):
    return log10(x)


# ----------------------------------------------------------------------
sec("0. PARAMETERS (Table 1)")
# ----------------------------------------------------------------------
p(f"Ne (diploid)              = {Ne:,.0f}")
p(f"gene copies 2Ne           = {Ncopies:,.0f} per generation")
p(f"mu per site per gen       = {mu}")
p(f"budget B                  = {B:,.0f} generations (6 Myr / 25 yr); generous {B_gen:,.0f} (9 Myr / 20 yr)")
p(f"log10(B)                  = {lg(B):.2f}")

# ----------------------------------------------------------------------
sec("A.1  CLASS I: SINGLE SUBSTITUTIONS  (Sec. 3.2, Appendix A.1)")
# ----------------------------------------------------------------------
# fixation probability: u = (1-exp(-2s))/(1-exp(-4 Ne s)) ~ 2s   (Kimura 1962)
s_fav = 0.01
u_fav = 2.0 * s_fav                                   # Haldane/Kimura small-s
t_fix_fav = (2.0 / s_fav) * math.log(2.0 * Ne)        # (2/s) ln(2Ne)
arrival_fav = 1.0 / (2.0 * Ne * mu * u_fav)           # 1/(2Ne mu u)
W1_fav = arrival_fav + t_fix_fav
p(f"favored allele s=0.01:  u = 2s = {u_fav}")
p(f"  arrival   1/(2Ne mu 2s) = {arrival_fav:,.3e} gen")
p(f"  fixation  (2/s) ln(2Ne) = {t_fix_fav:,.1f} gen")
p(f"  W1 total                  = {W1_fav:,.3e} gen  = {W1_fav*YR/1e6:.2f} Myr  ({W1_fav/B:.2f} B)")
p("")
neutral_specified_arrival = 1.0 / mu                  # 1/(2Ne mu * 1/(2Ne)) = 1/mu
t_fix_neu = 4.0 * Ne
p(f"neutral specified allele: arrival 1/mu = {neutral_specified_arrival:,.3e} gen; "
  f"conditional fixation 4Ne = {t_fix_neu:,.3e} gen ({t_fix_neu*YR/1e6:.1f} Myr)")
p("genome-wide neutral substitution rate = mu per site per generation (no arrival wait needed)")

# ----------------------------------------------------------------------
sec("A.3  SEGREGATING PROBABILITY OF A SPECIFIED k-TUPLE  (Eq. 4, Appendix A.3)")
# ----------------------------------------------------------------------
k = 8
C = 1.0e3
omega8 = C / 4.0 ** k
pseg = 2.0 * Ne * mu ** k * omega8 * (4.0 * Ne)
p(f"k=8, C_8 = 1e3  ->  omega_8 = C/4^8 = {omega8:.4e}")
p(f"P_seg(per locus) = 2Ne mu^8 omega_8 * 4Ne = {pseg:.3e}  (log10 = {lg(pseg):.2f})")
p(f"summed over 1e4 candidate loci: {pseg*1e4:.3e}  (log10 = {lg(pseg*1e4):.2f})")

# ----------------------------------------------------------------------
sec("3.3 / TABLE 3  FIRST-ARRIVAL TIME T = 4^k/(2Ne mu^k C_k)  (log10 generations)")
# ----------------------------------------------------------------------
p(f"{'Ne (diploid)':>14} {'2Ne':>10} | {'k=2 C=1':>9} {'k=5 C=1':>9} {'k=8 C=1':>9} | {'k=2 C=1e3':>10} {'k=5 C=1e3':>10} {'k=8 C=1e3':>10}")
for Ne_row, note in [(1e4, "hominin central"), (1e5, "ancestral hominoid peaks"), (1e6, "Lynch regime"),
                     (1e8, "unrealistic"), (1e10, "microbial")]:
    row = f"{Ne_row:>14,.0e} {2*Ne_row:>10,.0e} |"
    for C_ in (1.0, 1.0e3):
        for k_ in (2, 5, 8):
            row += f" {lg(T_arr(k_, C_, Ne_row)):>9.1f}"
        if C_ == 1.0:
            row += " |"
    p(row + f"   # {note}")
p()
p(f"budget log10 B = {lg(B):.2f}")
Ne_ext, C_ext = 1.0e9, 1.0e3
p(f"extreme case Ne=1e9, C_8=1e3, k=8: T = {T_arr(8, C_ext, Ne_ext):.3e} gen "
  f"(log10 = {lg(T_arr(8, C_ext, Ne_ext)):.2f}); excess over B = 10^{lg(T_arr(8, C_ext, Ne_ext)) - lg(B):.1f}")

# ----------------------------------------------------------------------
sec("PER-POSITION MULTIPLIER  (Sec. 3.3, Fig. 2 caption, Appendix B)")
# ----------------------------------------------------------------------
p(f"fixed C_k (omega = C/4^k): T(k+1)/T(k) = 4/mu      = {4.0/mu:.3e}   ({lg(4.0/mu):.2f} decades)")
p(f"fixed omega_k           : T(k+1)/T(k) = 1/mu      = {1.0/mu:.3e}   ({lg(1.0/mu):.2f} decades)")

# ----------------------------------------------------------------------
sec("SCENARIO 1  8-NUCLEOTIDE MICRO-RNA SEED MATCH, BOTH PRESITE REGIMES  (Sec. 6.2)")
# ----------------------------------------------------------------------
T_c1   = T_arr(8, 1.0)
T_c1e3 = T_arr(8, 1.0e3)
T_agg  = T_c1e3 / 1.0e10            # 1e3 motif variants x 1e4 loci x 1e3 windows
p("PRESITE-RICH (Durrett & Schmidt 2008, QUOTED):")
p("  a 1-kb region holds on average 3.85 six-of-eight partial matches to any target word;")
p("  a seven-of-eight presite exists with probability ~ 5/16.")
p("  8-mer completed from a 7-of-8 presite : 3.75e5 yr   (1.5e4 gen = 0.06 B)")
p("  PWM-flexible (any 7-of-8 binds)       : ~6e4 yr")
p("  de novo 6-mer in a 1-kb region        : ~1e5 yr")
p("  presite-free exact 8-mer (their sim.) : ~650 Myr  (= 108 B)")
p("")
p("PRESITE-FREE, single specified locus, COMPUTED via Eq. (3):")
p(f"  C_8 = 1   (omega_8 = 4^-8)      : T = {T_c1:.3e} gen  ({T_c1*YR:.2e} yr)")
p(f"  C_8 = 1e3 (omega_8 = 1.5e-2)    : T = {T_c1e3:.3e} gen  ({T_c1e3*YR:.2e} yr)")
p(f"  aggregate 1e10 target sites     : T = {T_agg:.3e} gen  -> excess over B = 10^{lg(T_agg)-lg(B):.1f}")

# ----------------------------------------------------------------------
sec("SCENARIO 2  DURRETT-SCHMIDT TWO-MUTATION REGULATORY CASE  (QUOTED, Table 2)")
# ----------------------------------------------------------------------
ds_mean_gen = 8.660258e6      # Theorem 1 prediction, Table 2 of DS 2008
ds_sim_gen  = 6.464920e6      # their simulation mean
p("parameters: n = 2e4 haploids; mu = 1e-8; u1 = 1e-7 (inactivation); u2 = 3.3e-9 (creation);")
p("            neutral intermediate; 1-kb region (DS 2008, Table 2)")
p(f"Theorem 1 mean = {ds_mean_gen:,.0f} generations = {ds_mean_gen*YR/1e6:.0f} Myr at 25 yr/gen")
p(f"simulation mean= {ds_sim_gen:,.0f} generations")
p(f"excess over B  = {ds_mean_gen/B:.1f} x")
p("")
p("CENSUS SENSITIVITY AXIS (arrival sums over census; establishment local to Ne; structured-population")
p("scaling linear in c = census/Ne; the panmictic alternative would set Ne = census, rejected by the")
p("three genealogical Ne measurements):")
for c in (1.0, 10.0, 100.0, 1000.0):
    Tc = ds_mean_gen / c
    p(f"  census {c:>5.0f} x Ne: mean = {Tc:,.3e} gen = {Tc*YR/1e6:.3f} Myr; "
      f"excess = {Tc/B:.4f} x B  ({'within budget' if Tc < B else 'outside budget'})")
p("plausibility anchor: habitable-area reconstructions for Pleistocene human populations place the")
p("Old World census at 2.1-8.3 million at the Last Glacial Maximum (Gautney & Holliday 2015,")
p("J Archaeol Sci 58:103-112), i.e., sustained census/Ne ratios of ~2e2-8e2 against Ne=1e4 - the")
p("hominin-specific 100x-1000x reading. Across other vertebrate taxa the sustained adult-census/Ne")
p("ratio is typically ~2-10 (Frankham 1995, Genet Res 66:95-107).")
p("establishment phase (fixation of the completed allele, ~4Ne = 1 Myr) is unaffected by census")
p("but is short relative to B in every reading.")
p("the presite-rich variants of the same biology (Scenario 1) run 1e4-6e4 gen and fit the budget")

# ----------------------------------------------------------------------
sec("A.2 / POLYNOMIAL TUNNELING CORRECTION AT THE DS PARAMETERS  (Sec. 3.5, Appendix A.2)")
# ----------------------------------------------------------------------
n_hap, u1_ds, u2_ds = 2.0e4, 1.0e-7, 3.3e-9
T_sim_ds = 1.0 / (n_hap * u1_ds * u2_ds)      # simultaneous arrival, both coordinated events
tunnel_factor = T_sim_ds / ds_mean_gen
p(f"simultaneous arrival for the DS two-event target: 1/(n u1 u2) = {T_sim_ds:.3e} gen")
p(f"DS tunneling (Theorem 1)                        = {ds_mean_gen:.3e} gen")
p(f"polynomial correction = {tunnel_factor:.3e}  (log10 = {lg(tunnel_factor):.2f}; ~10^4)")
p("worth ~10^4 at k=2 (decisive: brings the case to the budget margin), immaterial for k>=5")

# ----------------------------------------------------------------------
sec("SCENARIO 5  COMPOSITE OPERATION OF THE MECHANISM SET, PRESITE-FREE 8-MER  (Sec. 6.2)")
# ----------------------------------------------------------------------
T_raw   = T_c1e3            # single locus, C_8 = 1e3
T_sites = T_agg             # after 1e10 aggregate effective target sites
T_comp  = T_sites / tunnel_factor   # blanket tunneling credit at its DS (k=2) strength
p("best-case joint grants, presite-free 8-mer (C_8 = 1e3):")
p(f"  raw single locus                      : {T_raw:.3e} gen")
p(f"  (i) 1e10 aggregate effective sites    : {T_sites:.3e} gen")
p(f"  (ii) blanket tunneling x{tunnel_factor:.2e}      : {T_comp:.3e} gen")
p("  (iii) recombination/standing variation: inoperative (P_seg < 1e-52 over 1e4 loci; A.3)")
p("  (iv) IAD/genotype networks            : presuppose the seed (see K-S row, Table 2)")
p(f"composite ceiling  = {T_comp:.3e} gen = {T_comp*YR:.3e} yr")
p(f"composite excess over B = {T_comp/B:.3e}  (log10 = {lg(T_comp/B):.1f} -> ~10^{round(lg(T_comp/B))})")
p(f"compression achieved by joint operation = {T_raw/T_comp:.2e} (~{lg(T_raw/T_comp):.1f} orders)")
p("")
T_raw_k5_c1   = T_arr(5, 1.0)
T_raw_k5_c1e3 = T_arr(5, 1.0e3)
T_comp_k5_c1   = T_raw_k5_c1 / 1e10 / tunnel_factor
T_comp_k5_c1e3 = T_raw_k5_c1e3 / 1e10 / tunnel_factor
p("same joint grants applied at the class boundary k=5:")
p(f"  fully specified (C_5 = 1): raw excess = {T_raw_k5_c1/B:.2e} x B (log10 = {lg(T_raw_k5_c1/B):.1f});"
  f" composite excess = {T_comp_k5_c1/B:.2e} x B (log10 = {lg(T_comp_k5_c1/B):.1f})")
p(f"  maximal degeneracy (C_5 = 1e3): raw excess = {T_raw_k5_c1e3/B:.2e} x B (log10 = {lg(T_raw_k5_c1e3/B):.1f});"
  f" composite excess = {T_comp_k5_c1e3/B:.2e} x B (log10 = {lg(T_comp_k5_c1e3/B):.1f})")
p("class-level composite residual: >=1e18 x B at k=5 (fully specified), rising to ~1e41 x B at k=8")
p("Class IV composite route = Appendix A.4 supply model: marginally budget-compatible at")
p("generic-binding density (observed de novo gene rate), excluded at Axe-class specificity.")

# ----------------------------------------------------------------------
sec("SCENARIO 3  BEHE-SNOKE TWO-RESIDUE FEATURE VIA DUPLICATION  (Appendix A.5)")
# ----------------------------------------------------------------------
bs_N_req = 1.0e9              # their published threshold
bs_T_req = 1.0e8              # generations, their fixation window
bs_T_hominin = bs_T_req * (bs_N_req / Ne)   # linear extrapolation in 1/N, mutation-limited regime
p(f"published requirement: N >= {bs_N_req:.0e} to fix within {bs_T_req:.0e} generations")
p(f"linear extrapolation to Ne = {Ne:.0e}: T = {bs_T_hominin:.3e} gen = {bs_T_hominin*YR:.2e} yr")
p(f"excess over B = {bs_T_hominin/B:.2e} x")
p("Lynch (2005) rescue: ~1e6 yr 'with high probability' PROVIDED Ne > 1e6 (QUOTED; author notes")
p("the conclusion is 'a highly conservative conclusion with respect to both the time and")
p("population-size requirements'; his abstract: 'populations of only moderate size')")

# ----------------------------------------------------------------------
sec("SCENARIO 4 / A.4  DE NOVO BENEFICIAL ORF SUPPLY  (Sec. 5.6, Appendix A.4)")
# ----------------------------------------------------------------------
for Nc, n_tr, tag in [(1e5, 1e3, "central census"),
                     (1e6, 1e4, "generous census")]:
    trials = Nc * B * n_tr
    rho_fix = 60.0 / trials
    p(f"Nc={Nc:.0e}, n_tr={n_tr:.0e}: trials = {trials:.2e}; "
      f"required establishment-competent density rho_fix = {rho_fix:.2e}  [{tag}]")
    for s in (1e-2, 1e-3):
        est_cost = 1.0 / (2.0 * s)
        rho_orig = rho_fix * est_cost
        p(f"    s = {s:.0e}: establishment cost 1/(2s) = {est_cost:.0e}"
          f" -> required origination density = {rho_orig:.2e}")
p("benchmark densities (QUOTED):")
p("  Keefe-Szostak 2001: 4 ATP-binding 80-mers recovered from ~6e12 random sequences")
p("                      -> rho ~ 6.7e-13, quoted in the literature as order 1e-12 to 1e-11")
p("  Axe 2004           : rho <= 1e-77 for a specified 150-residue structured fold with function")
p("  -> origination requirement coincides with the generic-binding anchor at the generous census")
p("     reading (marginally supplyable), exceeds it at the central reading, and is tens of orders")
p("     ABOVE what Axe-class specificity would require (excluded if Axe-class) under every reading")
p("NOTE: the naive model T = 1/(2Ne mu^3 n_ORF rho) is NOT used in the manuscript; it double-counts")
p("specification by billing ORF origination itself to the standard model.")

# ----------------------------------------------------------------------
sec("TABLE 2  CLASS IV CONVERSION UNDER THE GENEROUS TRIAL CEILING  (theta = 1 novel chain/genome/gen)")
# ----------------------------------------------------------------------
theta = 1.0
for rho, label in [(1e-77, "Axe 2004 specified fold+function"),
                   (1e-11, "Keefe-Szostak any ATP binding (generous end)")]:
    trials_needed = 1.0 / rho
    T_gen = trials_needed / (2.0 * Ne * theta)
    p(f"rho = 1e{lg(rho):.0f}  [{label}]: trials = {trials_needed:.1e}; "
      f"T = {T_gen:.2e} gen = {T_gen*YR:.2e} yr; excess over B = 10^{lg(T_gen)-lg(B):.1f}")

# ----------------------------------------------------------------------
sec("SANFORD ET AL. 2015 (QUOTED)  MENDEL'S ACCOUNTANT, MODEL HOMININ POPULATION")
# ----------------------------------------------------------------------
p("parameters: N = 1e4; 20-yr generation; 50% selective elimination (their strong-selection")
p("setting - more extreme than empirical mortality schedules; treated as a sensitivity axis)")
p("shortest string tested k=1 : 1.53 Myr  (0.26 B)")
p("k=2                        : 84 Myr    (14 x B)")
p("k=5                        : 2.0 Gyr")
p("k=6                        : 4.24 Gyr")
p("k=8                        : 18.5 Gyr  (3.1e3 x B)")

# ----------------------------------------------------------------------
sec("FIGURE 2 DATA  T = 1/(2Ne mu^k), omega = 1  (waiting-time curves)")
# ----------------------------------------------------------------------
p("k=1 at Ne=1e4: T = 1/(2e4*1.25e-8) = 4.0e3 gen  (below B)")
for k_ in (2, 3, 5, 8):
    p(f"k={k_} at Ne=1e4: T = {T_arr(k_, 1.0):.2e} gen")
p("multiplier per added position (C fixed) = 4/mu = 3.2e8")

# ----------------------------------------------------------------------
sec("FIGURE 3 DATA  ESTIMATES vs BUDGET (years)")
# ----------------------------------------------------------------------
items = [
    ("Soft sweep from standing variant (conditional)",        5.0e4),
    ("8-mer from 7-of-8 presite (DS 2008, presite-rich)",     3.75e5),
    ("Neutral substitution, conditional fixation (4Ne)",      1.0e6),
    ("Specified favored substitution (arrival + fixation)",   W1_fav * YR),
    ("Two-mutation regulatory (DS), sustained 100x census",   ds_mean_gen / 100.0 * YR),
    ("Two-mutation regulatory (DS 2008 tunneling)",           216.0e6),
    ("Binding-only random protein (K-S, generous ceiling)",   5.0e6 * YR),
    ("2-residue feature via duplicate (B&S, extrapolated)",   bs_T_hominin * YR),
    ("Exact 8-nt string, presite-free, C=1e3 (this paper)",   T_c1e3 * YR),
    ("Composite best case, presite-free 8-mer (Scenario 5)",  T_comp * YR),
    ("Novel 150-aa fold at Axe density (generous ceiling)",   5.0e72 * YR),
]
for name, v in items:
    p(f"{v:>10.3e} yr   ({v/6.0e6:.2e} x budget)   {name}")

# ----------------------------------------------------------------------
sec("FIGURE 4 / MECHANISM MATRIX: no computation (literature matrix, Sec. 5)")
# ----------------------------------------------------------------------
p("Requirement matrix fig4_mechanism_matrix.png: qualitative literature synthesis; no numbers here.")

# ----------------------------------------------------------------------
# with open("/home/z/my-project/paper/numerical_outputs.txt", "w") as f:
#     f.write("\n".join(L) + "\n")
import pathlib
out_path = pathlib.Path(__file__).parent / "numerical_outputs.txt"
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print(f"\nsaved log -> {out_path}")
print("\nsaved log -> /home/z/my-project/paper/numerical_outputs.txt")


