module

import Mathlib.Analysis.Calculus.SmoothSeries
import Mathlib.Analysis.SpecialFunctions.SmoothTransition
public import Expdb.Mathlib.IteratedDeriv
import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
import Mathlib.Analysis.Calculus.ContDiff.FTaylorSeries
import Mathlib.Analysis.SpecificLimits.Basic

import Mathlib.Analysis.Calculus.ContDiff.Defs
import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
import Mathlib.Analysis.Calculus.TangentCone.Real
import Mathlib.Tactic.GCongr

/-!
# Smooth realization of jets and gluing on a closed interval

Adapted from `DifferentialGeometry/Analysis/Calculus/SmoothExtension/BorelHalfLine/Basic.lean`
and `DifferentialGeometry/Analysis/Calculus/SmoothExtension/JetGluing/Seam.lean` in
[`qinz1yang/differential-geometry`](https://github.com/qinz1yang/differential-geometry)
at commit `fd2f07fba4c6c102ba626c59153ad98da8f0d442`, licensed under Apache-2.0
(see the root LICENSE). Applicable upstream NOTICE:

> DifferentialGeometry
> Copyright 2026 The DifferentialGeometry contributors

Adaptations for Expdb: retained the Borel construction and closed-interval gluing theorem,
removed manifold dependencies, and added finite-order uniform bounds with a rescaling
argument for small prescribed derivatives. Construction helpers are private; namespaces
and proofs were updated for the current Mathlib API.
-/

@[expose] public section

namespace Expdb.SmoothExtension

noncomputable section
open Set Filter Topology
open scoped ContDiff

variable {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]

private def borelCutoff (s : ℝ) : ℝ := Real.smoothTransition (2 - s ^ 2)

private theorem borelCutoff_contDiff {n : ℕ∞} : ContDiff ℝ n borelCutoff := by
  unfold borelCutoff
  have h : ContDiff ℝ (n : WithTop ℕ∞) (fun s : ℝ => 2 - s ^ 2) :=
    contDiff_const.sub (contDiff_id.pow 2)
  exact Real.smoothTransition.contDiff.comp h

private theorem borelCutoff_eq_one {s : ℝ} (hs : s ^ 2 ≤ 1) : borelCutoff s = 1 := by
  refine Real.smoothTransition.one_of_one_le ?_
  linarith

private theorem borelCutoff_eq_zero {s : ℝ} (hs : 2 ≤ s ^ 2) : borelCutoff s = 0 := by
  refine Real.smoothTransition.zero_of_nonpos ?_
  linarith

private def borelBumpMono (n : ℕ) (s : ℝ) : ℝ := borelCutoff s * s ^ n / (Nat.factorial n : ℝ)

private theorem borelBumpMono_contDiff {N : ℕ∞} (n : ℕ) :
    ContDiff ℝ N (borelBumpMono n) :=
  ((borelCutoff_contDiff.mul ((contDiff_id).pow n)).div_const _)

private theorem borelBumpMono_comp_smul_contDiff {N : ℕ∞} (n : ℕ) (a : ℝ) :
    ContDiff ℝ N (fun s : ℝ => borelBumpMono n (a * s)) :=
  (borelBumpMono_contDiff n).comp (contDiff_const.mul contDiff_id)

private theorem borelBumpMono_support_subset (n : ℕ) :
    Function.support (borelBumpMono n) ⊆ Icc (-2 : ℝ) 2 := by
  intro s hs
  by_contra hmem
  apply hs
  have h2 : (2 : ℝ) ≤ s ^ 2 := by
    rcases not_and_or.1 hmem with h | h
    · have : s < -2 := lt_of_not_ge h
      nlinarith
    · have : (2 : ℝ) < s := lt_of_not_ge h
      nlinarith
  simp [borelBumpMono, borelCutoff_eq_zero h2]

private theorem borelBumpMono_hasCompactSupport (n : ℕ) :
    HasCompactSupport (borelBumpMono n) :=
  HasCompactSupport.of_support_subset_isCompact isCompact_Icc
    (borelBumpMono_support_subset n)

private theorem borelBumpMono_deriv_bound (n : ℕ) :
    ∃ G : ℝ, 0 ≤ G ∧ ∀ k ≤ n, ∀ s : ℝ,
      ‖iteratedDeriv k (borelBumpMono n) s‖ ≤ G := by
  have hbnd : ∀ k, ∃ C : ℝ, ∀ s : ℝ,
      ‖iteratedDeriv k (borelBumpMono n) s‖ ≤ C := by
    intro k
    have hcs : HasCompactSupport (iteratedFDeriv ℝ k (borelBumpMono n)) :=
      (borelBumpMono_hasCompactSupport n).iteratedFDeriv k
    have hcont : Continuous (iteratedFDeriv ℝ k (borelBumpMono n)) :=
      (borelBumpMono_contDiff (N := (⊤ : ℕ∞)) n).continuous_iteratedFDeriv
        (by exact_mod_cast le_top : (k : WithTop ℕ∞) ≤ ((⊤ : ℕ∞) : WithTop ℕ∞))
    obtain ⟨C, hC⟩ := hcs.exists_bound_of_continuous hcont
    refine ⟨C, fun s => ?_⟩
    rw [← norm_iteratedFDeriv_eq_norm_iteratedDeriv]
    exact hC s
  choose C hC using hbnd
  refine ⟨(Finset.range (n + 1)).sup' (by simp) (fun k => max (C k) 0), ?_, ?_⟩
  · exact le_trans (le_max_right _ _)
      (Finset.le_sup' (fun k => max (C k) 0) (by simp : (0 : ℕ) ∈ Finset.range (n + 1)))
  · intro k hk s
    refine le_trans (hC k s) ?_
    refine le_trans (le_max_left (C k) 0) ?_
    exact Finset.le_sup' (fun k => max (C k) 0)
      (Finset.mem_range.2 (Nat.lt_succ_of_le hk))


private def borelScale (c : ℕ → F) (n : ℕ) : ℝ :=
  max 1 (2 ^ n * Classical.choose (borelBumpMono_deriv_bound n) * (1 + ‖c n‖))

omit [NormedSpace ℝ F] in
private theorem one_le_borelScale (c : ℕ → F) (n : ℕ) : 1 ≤ borelScale c n :=
  le_max_left _ _

omit [NormedSpace ℝ F] in
private theorem borelScale_pos (c : ℕ → F) (n : ℕ) : 0 < borelScale c n :=
  lt_of_lt_of_le one_pos (one_le_borelScale c n)

private def borelTerm (c : ℕ → F) (n : ℕ) (t : ℝ) : F :=
  ((borelScale c n)⁻¹ ^ n * borelBumpMono n (borelScale c n * t)) • c n

private theorem borelTerm_contDiff {N : ℕ∞} (c : ℕ → F) (n : ℕ) :
    ContDiff ℝ N (borelTerm c n) := by
  refine ContDiff.smul ?_ contDiff_const
  exact contDiff_const.mul (borelBumpMono_comp_smul_contDiff n (borelScale c n))

private theorem borelTerm_hasCompactSupport (c : ℕ → F) (n : ℕ) :
    HasCompactSupport (borelTerm c n) := by
  apply HasCompactSupport.smul_right
    (f := fun t => (borelScale c n)⁻¹ ^ n * borelBumpMono n (borelScale c n * t))
  apply HasCompactSupport.of_support_subset_isCompact
    (K := Icc (-(2 / borelScale c n)) (2 / borelScale c n)) isCompact_Icc
  intro t ht
  by_contra hmem
  apply ht
  have hL0 : 0 < borelScale c n := borelScale_pos c n
  have h2 : 2 ≤ (borelScale c n * t) ^ 2 := by
    rcases not_and_or.1 hmem with h | h
    · have ht' : t < -(2 / borelScale c n) := lt_of_not_ge h
      have hlt : borelScale c n * t < borelScale c n * (-(2 / borelScale c n)) :=
        mul_lt_mul_of_pos_left ht' hL0
      rw [mul_neg, mul_div_cancel₀ _ hL0.ne'] at hlt
      nlinarith
    · have ht' : 2 / borelScale c n < t := lt_of_not_ge h
      have hlt : borelScale c n * (2 / borelScale c n) < borelScale c n * t :=
        mul_lt_mul_of_pos_left ht' hL0
      rw [mul_div_cancel₀ _ hL0.ne'] at hlt
      nlinarith
  simp only [borelBumpMono, borelCutoff_eq_zero h2, zero_mul, zero_div, mul_zero]

private theorem borelTerm_iteratedDeriv_global_bound (c : ℕ → F) (n k : ℕ) :
    ∃ B : ℝ, ∀ t : ℝ, ‖iteratedDeriv k (borelTerm c n) t‖ ≤ B := by
  have hcs : HasCompactSupport (iteratedFDeriv ℝ k (borelTerm c n)) :=
    (borelTerm_hasCompactSupport c n).iteratedFDeriv k
  have hcont : Continuous (iteratedFDeriv ℝ k (borelTerm c n)) :=
    (borelTerm_contDiff (N := (⊤ : ℕ∞)) c n).continuous_iteratedFDeriv
      (by exact_mod_cast le_top : (k : WithTop ℕ∞) ≤ ((⊤ : ℕ∞) : WithTop ℕ∞))
  obtain ⟨B, hB⟩ := hcs.exists_bound_of_continuous hcont
  refine ⟨B, fun t => ?_⟩
  rw [← norm_iteratedFDeriv_eq_norm_iteratedDeriv]
  exact hB t

private theorem borelTerm_iteratedDeriv (c : ℕ → F) (n k : ℕ) (t : ℝ) :
    iteratedDeriv k (borelTerm c n) t =
      ((borelScale c n)⁻¹ ^ n *
        (borelScale c n) ^ k * iteratedDeriv k (borelBumpMono n) (borelScale c n * t))
        • c n := by
  have hcd : ContDiffAt ℝ k (fun s => borelBumpMono n (borelScale c n * s)) t :=
    (borelBumpMono_comp_smul_contDiff (N := (k : ℕ∞)) n (borelScale c n)).contDiffAt
  have hsmul : iteratedDeriv k (borelTerm c n) t =
      (iteratedDeriv k (fun s => (borelScale c n)⁻¹ ^ n *
        borelBumpMono n (borelScale c n * s)) t) • c n :=
    iteratedDeriv_smul_const (contDiff_const.contDiffAt.mul hcd) _
  rw [hsmul]
  congr 1
  have hconst : iteratedDeriv k (fun s => (borelScale c n)⁻¹ ^ n *
      borelBumpMono n (borelScale c n * s)) t =
      (borelScale c n)⁻¹ ^ n •
        iteratedDeriv k (fun s => borelBumpMono n (borelScale c n * s)) t := by
    have := iteratedDeriv_const_smul
      (𝕜 := ℝ) (R := ℝ) (n := k) (x := t)
      (f := fun s => borelBumpMono n (borelScale c n * s)) hcd ((borelScale c n)⁻¹ ^ n)
    simpa only [Pi.smul_apply, smul_eq_mul] using! this
  rw [hconst]
  have hcomp : iteratedDeriv k (fun s => borelBumpMono n (borelScale c n * s)) t =
      (borelScale c n) ^ k • iteratedDeriv k (borelBumpMono n) (borelScale c n * t) := by
    have h := iteratedDeriv_comp_const_smul
      (n := k) (f := borelBumpMono n)
      (borelBumpMono_contDiff (N := (k : ℕ∞)) n) (borelScale c n)
    exact congrFun h t
  rw [hcomp]
  simp only [smul_eq_mul]
  ring

private theorem borelTerm_iteratedDeriv_bound (c : ℕ → F) {n k : ℕ} (hkn : k < n) (t : ℝ) :
    ‖iteratedDeriv k (borelTerm c n) t‖ ≤ (2 : ℝ)⁻¹ ^ n := by
  set G : ℝ := Classical.choose (borelBumpMono_deriv_bound n) with hGdef
  have hGspec := Classical.choose_spec (borelBumpMono_deriv_bound n)
  rw [← hGdef] at hGspec
  obtain ⟨hGnn, hGbound⟩ := hGspec
  set L : ℝ := borelScale c n with hLdef
  have hL1 : 1 ≤ L := one_le_borelScale c n
  have hL0 : 0 < L := borelScale_pos c n
  have hLmax : 2 ^ n * G * (1 + ‖c n‖) ≤ L := le_max_right _ _
  rw [borelTerm_iteratedDeriv, norm_smul]
  have hcn0 : (0 : ℝ) ≤ ‖c n‖ := norm_nonneg _
  have hD : ‖iteratedDeriv k (borelBumpMono n) (L * t)‖ ≤ G :=
    hGbound k (le_of_lt hkn) (L * t)
  have hpow : L⁻¹ ^ n * L ^ k ≤ L⁻¹ := by
    rw [inv_pow, inv_mul_le_iff₀ (by positivity)]
    have hrw : L ^ n * L⁻¹ = L ^ (n - 1) := by
      have : L ^ n = L ^ (n - 1) * L := by
        rw [← pow_succ, Nat.sub_add_cancel (by omega : 1 ≤ n)]
      rw [this, mul_assoc, mul_inv_cancel₀ hL0.ne', mul_one]
    rw [hrw]
    exact pow_le_pow_right₀ hL1 (by omega)
  have hcoeff : ‖L⁻¹ ^ n * L ^ k * iteratedDeriv k (borelBumpMono n) (L * t)‖
      ≤ L⁻¹ * G := by
    rw [norm_mul, norm_mul, norm_pow, norm_pow, Real.norm_eq_abs, Real.norm_eq_abs,
      abs_inv, abs_of_pos hL0]
    exact mul_le_mul hpow hD (norm_nonneg _) (by positivity)
  refine le_trans (mul_le_mul_of_nonneg_right hcoeff hcn0) ?_
  rcases eq_or_lt_of_le hGnn with hG0 | hGpos
  · rw [← hG0]
    simp only [mul_zero, zero_mul]
    positivity
  · have hLinv : L⁻¹ ≤ (2 ^ n * G * (1 + ‖c n‖))⁻¹ :=
      inv_anti₀ (by positivity) hLmax
    have h2pos : (0 : ℝ) < 2 ^ n := by positivity
    have h1cn : (0 : ℝ) < 1 + ‖c n‖ := by positivity
    have hcnle : ‖c n‖ ≤ 1 + ‖c n‖ := by linarith
    calc L⁻¹ * G * ‖c n‖
        ≤ (2 ^ n * G * (1 + ‖c n‖))⁻¹ * G * ‖c n‖ :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hLinv hGnn) hcn0
      _ = ‖c n‖ / (2 ^ n * (1 + ‖c n‖)) := by
          field_simp
      _ ≤ (2 : ℝ)⁻¹ ^ n := by
          rw [div_le_iff₀ (by positivity), inv_pow, inv_mul_eq_div, le_div_iff₀ h2pos]
          calc ‖c n‖ * 2 ^ n ≤ (1 + ‖c n‖) * 2 ^ n :=
                mul_le_mul_of_nonneg_right hcnle h2pos.le
            _ = 2 ^ n * (1 + ‖c n‖) := by ring

private def borelDomBound (c : ℕ → F) (k n : ℕ) : ℝ :=
  if k < n then (2 : ℝ)⁻¹ ^ n
  else Classical.choose (borelTerm_iteratedDeriv_global_bound c n k)

private theorem borelDomBound_summable (c : ℕ → F) (k : ℕ) :
    Summable (borelDomBound c k) := by
  have hgeom : Summable (fun n : ℕ => (2 : ℝ)⁻¹ ^ n) := by
    simpa only [one_div] using summable_geometric_two
  refine Summable.congr_cofinite hgeom ?_
  rw [Nat.cofinite_eq_atTop]
  filter_upwards [eventually_gt_atTop k] with n hn
  rw [borelDomBound, if_pos hn]

private theorem borelDomBound_bound (c : ℕ → F) (k n : ℕ) (t : ℝ) :
    ‖iteratedFDeriv ℝ k (borelTerm c n) t‖ ≤ borelDomBound c k n := by
  rw [norm_iteratedFDeriv_eq_norm_iteratedDeriv, borelDomBound]
  by_cases hkn : k < n
  · rw [if_pos hkn]
    exact borelTerm_iteratedDeriv_bound c hkn t
  · rw [if_neg hkn]
    exact Classical.choose_spec (borelTerm_iteratedDeriv_global_bound c n k) t

private theorem borelTerm_iteratedDeriv_zero (c : ℕ → F) (n k : ℕ) :
    iteratedDeriv k (borelTerm c n) 0 = (if n = k then (1 : ℝ) else 0) • c n := by
  have hL0 : 0 < borelScale c n := borelScale_pos c n
  have hev : borelTerm c n =ᶠ[𝓝 (0 : ℝ)]
      fun t => ((borelScale c n)⁻¹ ^ n * ((borelScale c n * t) ^ n / (Nat.factorial n : ℝ)))
        • c n := by
    refine eventuallyEq_of_mem (Ioo_mem_nhds (a := -(1 / borelScale c n))
      (b := 1 / borelScale c n) (by simpa using by positivity) (by positivity)) ?_
    intro t ht
    have h1 : (borelScale c n * t) ^ 2 ≤ 1 := by
      rw [mul_pow]
      have htabs : |t| < 1 / borelScale c n := by
        rw [abs_lt]; exact ⟨by simpa using ht.1, ht.2⟩
      have ht2 : t ^ 2 ≤ (1 / borelScale c n) ^ 2 := by
        rw [← sq_abs t]; exact pow_le_pow_left₀ (abs_nonneg t) htabs.le 2
      calc borelScale c n ^ 2 * t ^ 2
          ≤ borelScale c n ^ 2 * (1 / borelScale c n) ^ 2 :=
            mul_le_mul_of_nonneg_left ht2 (by positivity)
        _ = 1 := by field_simp
    simp only [borelTerm, borelBumpMono, borelCutoff_eq_one h1, one_mul]
  rw [Filter.EventuallyEq.iteratedDeriv_eq k hev]
  have hcd : ContDiffAt ℝ k
      (fun t => (borelScale c n)⁻¹ ^ n * ((borelScale c n * t) ^ n / (Nat.factorial n : ℝ))) 0 := by
    fun_prop
  rw [iteratedDeriv_smul_const hcd]
  have hpow : iteratedDeriv k
      (fun t => (borelScale c n)⁻¹ ^ n * ((borelScale c n * t) ^ n / (Nat.factorial n : ℝ))) 0 =
      if n = k then (1 : ℝ) else 0 := by
    have hfun : (fun t => (borelScale c n)⁻¹ ^ n *
        ((borelScale c n * t) ^ n / (Nat.factorial n : ℝ))) =
        (fun t => ((borelScale c n)⁻¹ ^ n * (borelScale c n) ^ n / (Nat.factorial n : ℝ)) *
          t ^ n) := by
      funext t; rw [mul_pow]; ring
    rw [hfun]
    simp only [iteratedDeriv_const_mul_field, iteratedDeriv_fun_pow_zero]
    have hLn : (borelScale c n)⁻¹ ^ n * borelScale c n ^ n = 1 := by
      rw [← mul_pow, inv_mul_cancel₀ hL0.ne', one_pow]
    rcases eq_or_ne n k with h | h
    · subst h
      rw [if_pos rfl, if_pos rfl]
      have hfac : (Nat.factorial n : ℝ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero n
      rw [div_mul_cancel₀ _ hfac, hLn]
    · rw [if_neg (fun h' => h h'.symm), if_neg h]
      simp
  rw [hpow]

private theorem borelTerm_bound_of_norm_le_one (n k : ℕ) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ (c : ℕ → F), ‖c n‖ ≤ 1 →
      ∀ t : ℝ, ‖iteratedDeriv k (borelTerm c n) t‖ ≤ B := by
  have hcs : HasCompactSupport (iteratedFDeriv ℝ k (borelBumpMono n)) :=
    (borelBumpMono_hasCompactSupport n).iteratedFDeriv k
  have hcont := (borelBumpMono_contDiff (N := (⊤ : ℕ∞)) n).continuous_iteratedFDeriv
    (by exact_mod_cast le_top : (k : WithTop ℕ∞) ≤ ((⊤ : ℕ∞) : WithTop ℕ∞))
  obtain ⟨D, hD⟩ := hcs.exists_bound_of_continuous hcont
  let G := Classical.choose (borelBumpMono_deriv_bound n)
  have hG : 0 ≤ G := (Classical.choose_spec (borelBumpMono_deriv_bound n)).1
  let U : ℝ := max 1 (2 ^ n * G * 2)
  have hU : 1 ≤ U := le_max_left _ _
  refine ⟨U ^ k * max D 0, by positivity, ?_⟩
  intro c hc t
  have hL : borelScale c n ≤ U := by
    apply max_le hU
    exact (mul_le_mul_of_nonneg_left (by linarith : 1 + ‖c n‖ ≤ 2)
      (by positivity : 0 ≤ 2 ^ n * G)).trans (le_max_right _ _)
  have hL0 := borelScale_pos c n
  have hi : (borelScale c n)⁻¹ ^ n ≤ 1 := by
    exact pow_le_one₀ (by positivity) (inv_le_one_of_one_le₀ (one_le_borelScale c n))
  have hd : ‖iteratedDeriv k (borelBumpMono n) (borelScale c n * t)‖ ≤ max D 0 := by
    rw [← norm_iteratedFDeriv_eq_norm_iteratedDeriv]
    exact (hD (borelScale c n * t)).trans (le_max_left _ _)
  rw [borelTerm_iteratedDeriv, norm_smul, norm_mul, norm_mul, norm_pow, norm_pow,
    Real.norm_of_nonneg (inv_nonneg.2 hL0.le), Real.norm_of_nonneg hL0.le]
  calc
    (borelScale c n)⁻¹ ^ n * borelScale c n ^ k *
        ‖iteratedDeriv k (borelBumpMono n) (borelScale c n * t)‖ * ‖c n‖
      ≤ 1 * U ^ k * max D 0 * 1 := by
        gcongr
    _ = _ := by ring

/-- Borel realization with bounds uniform in every jet beyond the controlled
finite prefix. This is needed when extending a family of model phases: high
endpoint derivatives may vary without a common bound. -/
private theorem exists_contDiff_iteratedDeriv_eq_with_bound [CompleteSpace F] (P : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (c : ℕ → F), (∀ n ≤ P, ‖c n‖ ≤ 1) →
      ∃ h : ℝ → F, ContDiff ℝ ∞ h ∧ (∀ k, iteratedDeriv k h 0 = c k) ∧
        ∀ k ≤ P, ∀ t, ‖iteratedDeriv k h t‖ ≤ C := by
  choose B hBnonneg hB using (borelTerm_bound_of_norm_le_one (F := F))
  let v (k n : ℕ) : ℝ := if n ≤ P then B n k else (2 : ℝ)⁻¹ ^ n
  have hvnonneg (k n : ℕ) : 0 ≤ v k n := by
    dsimp [v]
    split_ifs
    · exact hBnonneg n k
    · positivity
  have hvs (k : ℕ) : Summable (v k) := by
    refine Summable.congr_cofinite summable_geometric_two ?_
    rw [Nat.cofinite_eq_atTop]
    filter_upwards [eventually_gt_atTop P] with n hn
    simp [v, not_le.mpr hn]
  let C := ∑ k ∈ Finset.range (P + 1), ∑' n, v k n
  have hC : 0 ≤ C := Finset.sum_nonneg (fun k _ => tsum_nonneg (hvnonneg k))
  refine ⟨C, hC, ?_⟩
  intro c hc
  let h : ℝ → F := fun t => ∑' n, borelTerm c n t
  have hh : ContDiff ℝ ∞ h :=
    contDiff_tsum (N := (⊤ : ℕ∞)) (v := borelDomBound c)
      (fun n => borelTerm_contDiff c n)
      (fun k _ => borelDomBound_summable c k)
      (fun k n t _ => borelDomBound_bound c k n t)
  have hderiv (k : ℕ) (t : ℝ) :
      iteratedDeriv k h t = ∑' n, iteratedDeriv k (borelTerm c n) t := by
    have hFD := iteratedFDeriv_tsum_apply (N := (⊤ : ℕ∞)) (v := borelDomBound c)
      (fun n => borelTerm_contDiff c n)
      (fun k _ => borelDomBound_summable c k)
      (fun k n t _ => borelDomBound_bound c k n t)
      (le_top : (k : ℕ∞) ≤ (⊤ : ℕ∞)) t
    change iteratedDeriv k (fun t => ∑' n, borelTerm c n t) t = _
    rw [iteratedDeriv_eq_equiv_comp, Function.comp_apply, hFD,
      ← LinearIsometryEquiv.coe_toContinuousLinearEquiv,
      (ContinuousMultilinearMap.piFieldEquiv ℝ (Fin k) F).symm.toContinuousLinearEquiv.map_tsum]
    congr 1
  refine ⟨h, hh, ?_, ?_⟩
  · intro k
    rw [hderiv, tsum_congr (fun n => borelTerm_iteratedDeriv_zero c n k)]
    rw [tsum_eq_single k (fun n hn => by rw [if_neg hn, zero_smul])]
    rw [if_pos rfl, one_smul]
  · intro k hk t
    have hbound (n : ℕ) : ‖iteratedDeriv k (borelTerm c n) t‖ ≤ v k n := by
      dsimp [v]
      split_ifs with hn
      · exact hB n k c (hc n hn) t
      · exact borelTerm_iteratedDeriv_bound c (by omega) t
    have hnorms : Summable (fun n => ‖iteratedDeriv k (borelTerm c n) t‖) :=
      Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hbound (hvs k)
    rw [hderiv]
    exact ((norm_tsum_le_tsum_norm hnorms).trans
      (hnorms.tsum_le_tsum hbound (hvs k))).trans
        (Finset.single_le_sum (fun j _ => tsum_nonneg (hvnonneg j))
          (Finset.mem_range.2 (Nat.lt_succ_of_le hk)))

/-- Small finite jets admit an extension whose corresponding derivatives are
uniformly small. The constant depends only on the number of controlled jets. -/
theorem exists_contDiff_iteratedDeriv_eq_with_small_bound [CompleteSpace F] (P : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ δ : ℝ, 0 < δ → ∀ (c : ℕ → F),
      (∀ n ≤ P, ‖c n‖ ≤ δ) →
      ∃ h : ℝ → F, ContDiff ℝ ∞ h ∧ (∀ k, iteratedDeriv k h 0 = c k) ∧
        ∀ k ≤ P, ∀ t, ‖iteratedDeriv k h t‖ ≤ C * δ := by
  obtain ⟨C, hC, hrealize⟩ := exists_contDiff_iteratedDeriv_eq_with_bound (F := F) P
  refine ⟨C, hC, ?_⟩
  intro δ hδ c hc
  obtain ⟨h, hh, hjets, hbound⟩ := hrealize (fun n => δ⁻¹ • c n) (by
    intro n hn
    rw [norm_smul, Real.norm_of_nonneg (inv_nonneg.2 hδ.le)]
    calc
      δ⁻¹ * ‖c n‖ ≤ δ⁻¹ * δ := mul_le_mul_of_nonneg_left (hc n hn) (by positivity)
      _ = 1 := inv_mul_cancel₀ hδ.ne')
  refine ⟨fun t => δ • h t, contDiff_const.smul hh, ?_, ?_⟩
  · intro k
    rw [iteratedDeriv_fun_const_smul_field, hjets, smul_smul, mul_inv_cancel₀ hδ.ne',
      one_smul]
  · intro k hk t
    rw [iteratedDeriv_fun_const_smul_field, norm_smul, Real.norm_of_nonneg hδ.le]
    exact (mul_le_mul_of_nonneg_left (hbound k hk t) hδ.le).trans_eq (mul_comm _ _)

theorem contDiffOn_Icc_ite_of_jet_match {fL fR : ℝ → F} {a b x₀ : ℝ} (ha : a < x₀) (hb : x₀ < b)
    (hL : ContDiffOn ℝ ∞ fL (Icc a x₀)) (hR : ContDiffOn ℝ ∞ fR (Icc x₀ b))
    (hjet : ∀ n : ℕ, iteratedDerivWithin n fL (Icc a x₀) x₀
      = iteratedDerivWithin n fR (Icc x₀ b) x₀) :
    ContDiffOn ℝ ∞ (fun t => if t ≤ x₀ then fL t else fR t) (Icc a b) := by
  set f : ℝ → F := fun t => if t ≤ x₀ then fL t else fR t with hf_def
  set pL : ℝ → FormalMultilinearSeries ℝ ℝ F := ftaylorSeriesWithin ℝ fL (Icc a x₀) with hpL_def
  set pR : ℝ → FormalMultilinearSeries ℝ ℝ F := ftaylorSeriesWithin ℝ fR (Icc x₀ b) with hpR_def
  set p : ℝ → FormalMultilinearSeries ℝ ℝ F :=
    fun t => if t ≤ x₀ then pL t else pR t with hp_def
  have hTL : HasFTaylorSeriesUpToOn ∞ fL pL (Icc a x₀) :=
    hL.ftaylorSeriesWithin (uniqueDiffOn_Icc ha)
  have hTR : HasFTaylorSeriesUpToOn ∞ fR pR (Icc x₀ b) :=
    hR.ftaylorSeriesWithin (uniqueDiffOn_Icc hb)
  have hjetF : ∀ n : ℕ, pL x₀ n = pR x₀ n := by
    intro n
    have heL := congrFun (iteratedFDerivWithin_eq_equiv_comp (𝕜 := ℝ) (n := n) (f := fL)
      (s := Icc a x₀)) x₀
    have heR := congrFun (iteratedFDerivWithin_eq_equiv_comp (𝕜 := ℝ) (n := n) (f := fR)
      (s := Icc x₀ b)) x₀
    simp only [Function.comp_apply] at heL heR
    have hLR : iteratedFDerivWithin ℝ n fL (Icc a x₀) x₀
        = iteratedFDerivWithin ℝ n fR (Icc x₀ b) x₀ := by
      rw [heL, heR, hjet n]
    simpa only [hpL_def, hpR_def, ftaylorSeriesWithin] using hLR
  have hEqL : ∀ m : ℕ, Set.EqOn (fun y => p y m) (fun y => pL y m) (Icc a x₀) := by
    intro m y hy
    simp only [hp_def, if_pos hy.2]
  have hEqR : ∀ m : ℕ, Set.EqOn (fun y => p y m) (fun y => pR y m) (Icc x₀ b) := by
    intro m y hy
    rcases eq_or_lt_of_le hy.1 with hy0 | hy0
    · subst hy0
      simp only [hp_def, if_pos (le_refl x₀), hjetF]
    · simp only [hp_def, if_neg (not_le.mpr hy0)]
  have hzero : ∀ x ∈ Icc a b, (p x 0).curry0 = f x := by
    intro x hx
    by_cases hx' : x ≤ x₀
    · have hval : (pL x 0).curry0 = fL x := hTL.zero_eq x ⟨hx.1, hx'⟩
      simp only [hp_def, hf_def, if_pos hx', hval]
    · have hx'' : x₀ ≤ x := le_of_lt (not_le.mp hx')
      have hval : (pR x 0).curry0 = fR x := hTR.zero_eq x ⟨hx'', hx.2⟩
      simp only [hp_def, hf_def, if_neg hx', hval]
  have hm_lt : ∀ m : ℕ, (m : WithTop ℕ∞) < ∞ := fun m => by
    exact_mod_cast (Nat.cast_lt.mpr m.lt_succ_self).trans_le le_top
  have hderiv : ∀ m : ℕ, ∀ x ∈ Icc a b,
      HasFDerivWithinAt (fun y => p y m) (p x m.succ).curryLeft (Icc a b) x := by
    intro m x hx
    have hpL_succ : ∀ y : ℝ, y ≤ x₀ → p y m.succ = pL y m.succ := fun y hy => by
      simp only [hp_def, if_pos hy]
    have hpR_succ : ∀ y : ℝ, x₀ ≤ y → p y m.succ = pR y m.succ := by
      intro y hy
      rcases eq_or_lt_of_le hy with hy0 | hy0
      · subst hy0; simp only [hp_def, if_pos (le_refl x₀), hjetF]
      · simp only [hp_def, if_neg (not_le.mpr hy0)]
    rcases lt_trichotomy x x₀ with hlt | heq | hgt
    · have hdL : HasFDerivWithinAt (fun y => pL y m) (pL x m.succ).curryLeft (Icc a x₀) x :=
        hTL.fderivWithin m (hm_lt m) x ⟨hx.1, hlt.le⟩
      have hmem : Icc a x₀ ∈ 𝓝[Icc a b] x :=
        mem_of_superset
          (inter_mem (nhdsWithin_le_nhds (isOpen_Iio.mem_nhds hlt)) self_mem_nhdsWithin)
          (fun y hy => ⟨hy.2.1, hy.1.le⟩)
      have hdL' : HasFDerivWithinAt (fun y => pL y m) (pL x m.succ).curryLeft (Icc a b) x :=
        hdL.mono_of_mem_nhdsWithin hmem
      have hee : (fun y => p y m) =ᶠ[𝓝[Icc a b] x] (fun y => pL y m) :=
        eventuallyEq_iff_exists_mem.mpr ⟨Icc a x₀, hmem, fun y hy => hEqL m hy⟩
      have hfd : HasFDerivWithinAt (fun y => p y m) (pL x m.succ).curryLeft (Icc a b) x :=
        hdL'.congr_of_eventuallyEq hee (hEqL m ⟨hx.1, hlt.le⟩)
      rw [hpL_succ x hlt.le]
      exact hfd
    · subst x
      have hdL0 : HasFDerivWithinAt (fun y => pL y m) (pL x₀ m.succ).curryLeft (Icc a x₀) x₀ :=
        hTL.fderivWithin m (hm_lt m) x₀ ⟨ha.le, le_refl x₀⟩
      have hdL0' : HasFDerivWithinAt (fun y => p y m) (pL x₀ m.succ).curryLeft (Icc a x₀) x₀ :=
        hdL0.congr (hEqL m) (hEqL m ⟨ha.le, le_refl x₀⟩)
      have hdR0 : HasFDerivWithinAt (fun y => pR y m) (pR x₀ m.succ).curryLeft (Icc x₀ b) x₀ :=
        hTR.fderivWithin m (hm_lt m) x₀ ⟨le_refl x₀, hb.le⟩
      have hdR0' : HasFDerivWithinAt (fun y => p y m) (pL x₀ m.succ).curryLeft (Icc x₀ b) x₀ := by
        have hval : (pR x₀ m.succ).curryLeft = (pL x₀ m.succ).curryLeft := by rw [hjetF]
        rw [← hval]
        exact hdR0.congr (hEqR m) (hEqR m ⟨le_refl x₀, hb.le⟩)
      have hunion : HasFDerivWithinAt (fun y => p y m) (pL x₀ m.succ).curryLeft
          (Icc a x₀ ∪ Icc x₀ b) x₀ := hdL0'.union hdR0'
      rw [Set.Icc_union_Icc_eq_Icc ha.le hb.le] at hunion
      rw [hpL_succ x₀ (le_refl x₀)]
      exact hunion
    · have hdR : HasFDerivWithinAt (fun y => pR y m) (pR x m.succ).curryLeft (Icc x₀ b) x :=
        hTR.fderivWithin m (hm_lt m) x ⟨hgt.le, hx.2⟩
      have hmem : Icc x₀ b ∈ 𝓝[Icc a b] x :=
        mem_of_superset
          (inter_mem (nhdsWithin_le_nhds (isOpen_Ioi.mem_nhds hgt)) self_mem_nhdsWithin)
          (fun y hy => ⟨hy.1.le, hy.2.2⟩)
      have hdR' : HasFDerivWithinAt (fun y => pR y m) (pR x m.succ).curryLeft (Icc a b) x :=
        hdR.mono_of_mem_nhdsWithin hmem
      have hee : (fun y => p y m) =ᶠ[𝓝[Icc a b] x] (fun y => pR y m) :=
        eventuallyEq_iff_exists_mem.mpr ⟨Icc x₀ b, hmem, fun y hy => hEqR m hy⟩
      have hfd : HasFDerivWithinAt (fun y => p y m) (pR x m.succ).curryLeft (Icc a b) x :=
        hdR'.congr_of_eventuallyEq hee (hEqR m ⟨hgt.le, hx.2⟩)
      rw [hpR_succ x hgt.le]
      exact hfd
  have hTaylor : HasFTaylorSeriesUpToOn ∞ f p (Icc a b) :=
    (hasFTaylorSeriesUpToOn_top_iff' (le_refl _)).mpr ⟨hzero, hderiv⟩
  exact hTaylor.contDiffOn

/-! ## Controlled extension on a closed interval

The following extension theorem is proved within Expdb. -/

/-- A smooth function on `[1, 2]` extends exactly to any larger compact interval,
with finite derivative bounds uniform over the original function. -/
theorem exists_contDiffOn_Icc_extension_with_bound (P : ℕ) :
    ∃ C : ℝ, 1 ≤ C ∧ ∀ δ : ℝ, 0 < δ → ∀ f : ℝ → ℝ,
      ContDiffOn ℝ ∞ f (Icc 1 2) →
      (∀ k ≤ P, ∀ x ∈ Icc 1 2, ‖iteratedDerivWithin k f (Icc 1 2) x‖ ≤ δ) →
      ∀ l u : ℝ, l < 1 → 2 < u →
      ∃ g : ℝ → ℝ, ContDiffOn ℝ ∞ g (Icc l u) ∧ EqOn g f (Icc 1 2) ∧
        ∀ k ≤ P, ∀ x ∈ Icc l u, ‖iteratedDerivWithin k g (Icc l u) x‖ ≤ C * δ := by
  obtain ⟨C, hC, hrealize⟩ :=
    exists_contDiff_iteratedDeriv_eq_with_small_bound (F := ℝ) P
  refine ⟨max 1 C, le_max_left _ _, ?_⟩
  intro δ hδ f hf hbound l u hl hu
  obtain ⟨L, hL, hLj, hLb⟩ := hrealize δ hδ
    (fun k => iteratedDerivWithin k f (Icc 1 2) 1)
    (fun k hk => hbound k hk 1 (by norm_num))
  obtain ⟨R, hR, hRj, hRb⟩ := hrealize δ hδ
    (fun k => iteratedDerivWithin k f (Icc 1 2) 2)
    (fun k hk => hbound k hk 2 (by norm_num))
  let fL : ℝ → ℝ := fun x => L (x - 1)
  let fR : ℝ → ℝ := fun x => R (x - 2)
  have hfL : ContDiff ℝ ∞ fL := hL.comp (contDiff_id.sub contDiff_const)
  have hfR : ContDiff ℝ ∞ fR := hR.comp (contDiff_id.sub contDiff_const)
  have hLderiv (k : ℕ) (x : ℝ) : iteratedDeriv k fL x = iteratedDeriv k L (x - 1) :=
    congrFun (iteratedDeriv_comp_sub_const k L 1) x
  have hRderiv (k : ℕ) (x : ℝ) : iteratedDeriv k fR x = iteratedDeriv k R (x - 2) :=
    congrFun (iteratedDeriv_comp_sub_const k R 2) x
  have hLzero : fL 1 = f 1 := by
    simpa [fL, iteratedDeriv_zero, iteratedDerivWithin_zero] using hLj 0
  have hRzero : fR 2 = f 2 := by
    simpa [fR, iteratedDeriv_zero, iteratedDerivWithin_zero] using hRj 0
  let g₁ : ℝ → ℝ := fun x => if x ≤ 1 then fL x else f x
  have hg₁ : ContDiffOn ℝ ∞ g₁ (Icc l 2) := by
    apply contDiffOn_Icc_ite_of_jet_match hl (by norm_num)
      hfL.contDiffOn hf
    intro k
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hl)
      (hfL.of_le (by exact_mod_cast (le_top : (k : ℕ∞) ≤ ⊤))).contDiffAt
        ⟨hl.le, le_rfl⟩, hLderiv, sub_self, hLj]
  have hg₁eq : EqOn g₁ f (Icc 1 2) := by
    intro x hx
    dsimp [g₁]
    split_ifs with h
    · have : x = 1 := le_antisymm h hx.1
      subst x
      exact hLzero
    · rfl
  have hsub : Icc 1 2 ⊆ Icc l 2 := fun x hx => ⟨hl.le.trans hx.1, hx.2⟩
  let g : ℝ → ℝ := fun x => if x ≤ 2 then g₁ x else fR x
  have hg : ContDiffOn ℝ ∞ g (Icc l u) := by
    apply contDiffOn_Icc_ite_of_jet_match (by linarith) hu
      hg₁ hfR.contDiffOn
    intro k
    rw [iteratedDerivWithin_eq_of_restriction hsub (uniqueDiffOn_Icc (by norm_num))
      (uniqueDiffOn_Icc (by linarith)) hg₁ hg₁eq (by norm_num) k]
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hu)
      (hfR.of_le (by exact_mod_cast (le_top : (k : ℕ∞) ≤ ⊤))).contDiffAt
        ⟨le_rfl, hu.le⟩, hRderiv, sub_self, hRj]
  have hgeq : EqOn g f (Icc 1 2) := by
    intro x hx
    exact (if_pos hx.2).trans (hg₁eq hx)
  have hgLeq : EqOn g fL (Icc l 1) := by
    intro x hx
    dsimp [g, g₁]
    rw [if_pos (hx.2.trans (by norm_num)), if_pos hx.2]
  have hgReq : EqOn g fR (Icc 2 u) := by
    intro x hx
    dsimp [g]
    split_ifs with h
    · have : x = 2 := le_antisymm h hx.1
      subst x
      exact (hg₁eq (by norm_num)).trans hRzero.symm
    · rfl
  have hJ : UniqueDiffOn ℝ (Icc l u) := uniqueDiffOn_Icc (by linarith)
  have hleftSub : Icc l 1 ⊆ Icc l u := fun y hy => ⟨hy.1, by linarith [hy.2]⟩
  have hrightSub : Icc 2 u ⊆ Icc l u := fun y hy => ⟨by linarith [hy.1], hy.2⟩
  have hmiddleSub : Icc 1 2 ⊆ Icc l u :=
    fun y hy => ⟨hl.le.trans hy.1, hy.2.trans hu.le⟩
  refine ⟨g, hg, hgeq, ?_⟩
  intro k hk x hx
  by_cases hleft : x ≤ 1
  · have hxL : x ∈ Icc l 1 := ⟨hx.1, hleft⟩
    rw [iteratedDerivWithin_eq_of_restriction
      hleftSub
      (uniqueDiffOn_Icc hl) hJ hg hgLeq hxL k,
      iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hl)
        (hfL.of_le (by exact_mod_cast (le_top : (k : ℕ∞) ≤ ⊤))).contDiffAt hxL, hLderiv]
    exact (hLb k hk _).trans (mul_le_mul_of_nonneg_right (le_max_right _ _) hδ.le)
  · by_cases hright : 2 ≤ x
    · have hxR : x ∈ Icc 2 u := ⟨hright, hx.2⟩
      rw [iteratedDerivWithin_eq_of_restriction
        hrightSub
        (uniqueDiffOn_Icc hu) hJ hg hgReq hxR k,
        iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hu)
          (hfR.of_le (by exact_mod_cast (le_top : (k : ℕ∞) ≤ ⊤))).contDiffAt hxR, hRderiv]
      exact (hRb k hk _).trans (mul_le_mul_of_nonneg_right (le_max_right _ _) hδ.le)
    · have hxF : x ∈ Icc 1 2 := ⟨(lt_of_not_ge hleft).le, (lt_of_not_ge hright).le⟩
      rw [iteratedDerivWithin_eq_of_restriction
        hmiddleSub
        (uniqueDiffOn_Icc (by norm_num)) hJ hg hgeq hxF k]
      exact (hbound k hk x hxF).trans (by
        simpa using mul_le_mul_of_nonneg_right (le_max_left (1 : ℝ) C) hδ.le)

/-! ## Smooth multiplication by a cutoff

The following lemma is proved within Expdb. -/

/-- A smooth cutoff times a function smooth on an open neighbourhood of its support is
smooth. -/
theorem contDiff_mul_of_tsupport {χ g : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U)
    (hχ : ContDiff ℝ ∞ χ) (hg : ContDiffOn ℝ ∞ g U) (hsupp : tsupport χ ⊆ U) :
    ContDiff ℝ ∞ (fun x ↦ χ x * g x) := by
  rw [contDiff_iff_contDiffAt]
  intro x
  by_cases hx : x ∈ U
  · exact hχ.contDiffAt.mul (hg.contDiffAt (hU.mem_nhds hx))
  · have hnot : x ∉ tsupport χ := fun h ↦ hx (hsupp h)
    have hzero : (fun y ↦ χ y * g y) =ᶠ[𝓝 x] fun _ ↦ 0 := by
      filter_upwards [(notMem_tsupport_iff_eventuallyEq.1 hnot)] with y hy
      simp [show χ y = 0 from hy]
    exact (contDiffAt_const (c := (0 : ℝ))).congr_of_eventuallyEq hzero


/-- A function smooth on a larger interval has a globally smooth representative on the
smaller positive interval. -/
theorem exists_contDiff_eqOn_Icc (l u : ℝ) (hl : 0 < l) (hlu : l < u)
    (g : ℝ → ℝ) (hg : ContDiffOn ℝ ∞ g (Icc (l / 2) (2 * u))) :
    ∃ G : ℝ → ℝ, ContDiff ℝ ∞ G ∧ EqOn G g (Icc l u) := by
  let χ : ContDiffBump ((l + u) / 2 : ℝ) :=
    ⟨(u - l) / 2 + l / 8, (u - l) / 2 + l / 4, by linarith, by linarith⟩
  let G : ℝ → ℝ := fun x => χ x * g x
  have hsupport : tsupport (χ : ℝ → ℝ) ⊆ Ioo (l / 2) (2 * u) := by
    intro x hx
    rw [χ.tsupport_eq, Metric.mem_closedBall, Real.dist_eq, abs_le] at hx
    change -((u - l) / 2 + l / 4) ≤ x - (l + u) / 2 ∧
      x - (l + u) / 2 ≤ (u - l) / 2 + l / 4 at hx
    constructor <;> linarith [hx.1, hx.2]
  have hG : ContDiff ℝ ∞ G :=
    contDiff_mul_of_tsupport isOpen_Ioo χ.contDiff
      (hg.mono Ioo_subset_Icc_self) hsupport
  refine ⟨G, hG, ?_⟩
  intro x hx
  have hχ : χ x = 1 := by
    apply χ.one_of_mem_closedBall
    rw [Metric.mem_closedBall, Real.dist_eq, abs_le]
    change -((u - l) / 2 + l / 8) ≤ x - (l + u) / 2 ∧
      x - (l + u) / 2 ≤ (u - l) / 2 + l / 8
    constructor <;> linarith [hx.1, hx.2]
  simp [G, hχ]

end
end Expdb.SmoothExtension
