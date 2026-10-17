module

public import Expdb.ExponentialSums.PhaseFunctions
public import Expdb.Mathlib.SmoothExtension
public import Expdb.Mathlib.IteratedDeriv
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Extension of approximate model phases

This file extends approximate model phases to larger positive intervals while
preserving their values on `[1, 2]` and controlling finitely many derivative errors.
-/

@[expose] public section

namespace Expdb

open Set
open scoped ContDiff Topology
noncomputable section

/-- A primitive of the reference derivative on the positive half-line. -/
private def modelPhasePrimitive (σ : ℝ) (x : ℝ) : ℝ :=
  if σ = 1 then Real.log x else x ^ (1 - σ) / (1 - σ)

private theorem hasDerivAt_modelPhasePrimitive (σ : ℝ) {x : ℝ} (hx : 0 < x) :
    HasDerivAt (modelPhasePrimitive σ) (modelPhase σ x) x := by
  change HasDerivAt (fun y => if σ = 1 then Real.log y else y ^ (1 - σ) / (1 - σ))
    (modelPhase σ x) x
  by_cases hσ : σ = 1
  · subst σ
    simpa [modelPhasePrimitive, modelPhase, Real.rpow_neg_one] using Real.hasDerivAt_log hx.ne'
  · have h := (Real.hasDerivAt_rpow_const (p := 1 - σ) (Or.inl hx.ne')).div_const (1 - σ)
    have hd : (1 - σ) * x ^ (1 - σ - 1) / (1 - σ) = x ^ (-σ) := by
      have hs : 1 - σ ≠ 0 := sub_ne_zero.mpr (Ne.symm hσ)
      rw [show 1 - σ - 1 = -σ by ring, mul_div_cancel_left₀ _ hs]
    simpa only [modelPhasePrimitive, if_neg hσ, modelPhase, hd] using h

private theorem contDiffOn_modelPhasePrimitive (σ : ℝ) {s : Set ℝ}
    (hs : ∀ x ∈ s, 0 < x) : ContDiffOn ℝ ∞ (modelPhasePrimitive σ) s := by
  intro x hx
  change ContDiffWithinAt ℝ ∞
    (fun y => if σ = 1 then Real.log y else y ^ (1 - σ) / (1 - σ)) s x
  by_cases hσ : σ = 1
  · simp only [hσ]
    exact (Real.contDiffAt_log.mpr (hs x hx).ne').contDiffWithinAt
  · simp only [if_neg hσ]
    exact ((Real.contDiffAt_rpow_const_of_ne (hs x hx).ne').div_const _).contDiffWithinAt

private theorem iteratedDerivWithin_modelPhasePrimitive (σ c : ℝ) {s : Set ℝ}
    (hs : UniqueDiffOn ℝ s) (hpos : ∀ x ∈ s, 0 < x) (p : ℕ) {x : ℝ} (hx : x ∈ s) :
    iteratedDerivWithin (p + 1) (fun y ↦ c + modelPhasePrimitive σ y) s x =
      iteratedDerivWithin p (modelPhase σ) s x := by
  rw [iteratedDerivWithin_const_add (by omega), iteratedDerivWithin_succ']
  apply iteratedDerivWithin_congr (s := s) _ hx
  intro y hy
  exact (hasDerivAt_modelPhasePrimitive σ (hpos y hy)).hasDerivWithinAt.derivWithin
    (hs.uniqueDiffWithinAt hy)

/-- An approximate model phase extends exactly to any larger positive compact
interval. The model derivative errors remain uniformly small through the same
order; the constant is independent of the phase and of its higher derivatives. -/
theorem exists_approximateModelPhase_extension (P : ℕ) :
    ∃ C : ℝ, 1 ≤ C ∧ ∀ δ : ℝ, 0 < δ → ∀ (σ : ℝ) (F : ℝ → ℝ),
      IsApproximateModelPhaseFunction F σ P δ →
      ∀ l u : ℝ, 0 < l → l < 1 → 2 < u →
      ∃ G : ℝ → ℝ, ContDiffOn ℝ ∞ G (Icc l u) ∧ EqOn G F phaseInterval ∧
        ∀ p ≤ P, ∀ x ∈ Icc l u,
          ‖iteratedDerivWithin (p + 1) G (Icc l u) x -
            iteratedDerivWithin p (modelPhase σ) (Icc l u) x‖ ≤ C * δ := by
  obtain ⟨C, hC, hextend⟩ := SmoothExtension.exists_contDiffOn_Icc_extension_with_bound (P + 1)
  refine ⟨C, hC, ?_⟩
  intro δ hδ σ F hF l u hl0 hl hu
  have hpos : ∀ x ∈ phaseInterval, 0 < x := fun x hx => zero_lt_one.trans_le hx.1
  have hQ := contDiffOn_modelPhasePrimitive σ hpos
  let c : ℝ := F 1 - modelPhasePrimitive σ 1
  let Q : ℝ → ℝ := fun x => c + modelPhasePrimitive σ x
  let f : ℝ → ℝ := fun x => F x - Q x
  have hQc : ContDiffOn ℝ ∞ Q phaseInterval := contDiffOn_const.add hQ
  have hf : ContDiffOn ℝ ∞ f phaseInterval := hF.1.sub hQc
  have hf0 : f 1 = 0 := by dsimp [f, Q, c]; ring
  have herror (p : ℕ) {x : ℝ} (hx : x ∈ phaseInterval) :
      iteratedDerivWithin (p + 1) f phaseInterval x = modelPhaseErrorAt F σ p x := by
    change iteratedDerivWithin (p + 1) (F - Q) phaseInterval x = _
    rw [iteratedDerivWithin_sub hx uniqueDiffOn_phaseInterval
      ((hF.1.of_le (by exact_mod_cast (le_top : ((p + 1) : ℕ∞) ≤ ⊤))) x hx)
      ((hQc.of_le (by exact_mod_cast (le_top : ((p + 1) : ℕ∞) ≤ ⊤))) x hx)]
    change iteratedDerivWithin (p + 1) F phaseInterval x -
      iteratedDerivWithin (p + 1) (fun x => c + modelPhasePrimitive σ x) phaseInterval x = _
    rw [iteratedDerivWithin_modelPhasePrimitive σ c uniqueDiffOn_phaseInterval hpos p hx]
    rfl
  have hfirst : ∀ x ∈ phaseInterval, ‖derivWithin f phaseInterval x‖ ≤ δ := by
    intro x hx
    have h := hF.2 0 (by omega) ⟨x, hx⟩
    rw [← herror 0 hx] at h
    simpa only [Nat.zero_add, iteratedDerivWithin_one] using h
  have hb : ∀ k ≤ P + 1, ∀ x ∈ phaseInterval,
      ‖iteratedDerivWithin k f phaseInterval x‖ ≤ δ := by
    intro k hk x hx
    cases k with
    | zero =>
      have h := (convex_Icc (1 : ℝ) 2).norm_image_sub_le_of_norm_derivWithin_le
        (hf.differentiableOn (by simp)) hfirst
        (show (1 : ℝ) ∈ phaseInterval by norm_num [phaseInterval]) hx
      have hnorm : ‖x - 1‖ ≤ 1 := by
        rw [Real.norm_of_nonneg (sub_nonneg.2 hx.1)]
        linarith [hx.2]
      simp only [hf0, sub_zero] at h
      exact h.trans (by simpa using mul_le_mul_of_nonneg_left hnorm hδ.le)
    | succ p =>
      rw [herror p hx]
      exact hF.2 p (by omega) ⟨x, hx⟩
  obtain ⟨g, hg, hgeq, hgb⟩ := hextend δ hδ f hf hb l u hl hu
  let G : ℝ → ℝ := fun x => g x + Q x
  have hposJ : ∀ x ∈ Icc l u, 0 < x := fun x hx => hl0.trans_le hx.1
  have hQJ : ContDiffOn ℝ ∞ Q (Icc l u) :=
    contDiffOn_const.add (contDiffOn_modelPhasePrimitive σ hposJ)
  refine ⟨G, hg.add hQJ, ?_, ?_⟩
  · intro x hx
    dsimp [G]
    rw [hgeq hx]
    dsimp [f]
    ring
  · intro p hp x hx
    have hJ := uniqueDiffOn_Icc (show l < u by linarith)
    rw [iteratedDerivWithin_fun_add hx hJ
      ((hg.of_le (by exact_mod_cast (le_top : ((p + 1) : ℕ∞) ≤ ⊤))) x hx)
      ((hQJ.of_le (by exact_mod_cast (le_top : ((p + 1) : ℕ∞) ≤ ⊤))) x hx)]
    change ‖iteratedDerivWithin (p + 1) g (Icc l u) x +
      iteratedDerivWithin (p + 1) (fun x => c + modelPhasePrimitive σ x) (Icc l u) x - _‖ ≤ _
    rw [iteratedDerivWithin_modelPhasePrimitive σ c hJ hposJ p hx, add_sub_cancel_right]
    exact hgb (p + 1) (by omega) x hx

/-- The controlled extension can be chosen globally smooth. This allows the
ordinary-derivative normal-form machinery to be applied on a larger positive
interval, without imposing any conditions on the original phase outside `[1,2]`. -/
theorem exists_global_approximateModelPhase_extension (P : ℕ) :
    ∃ C : ℝ, 1 ≤ C ∧ ∀ δ : ℝ, 0 < δ → ∀ σ : ℝ, ∀ F : ℝ → ℝ,
      IsApproximateModelPhaseFunction F σ P δ →
      ∀ l u : ℝ, 0 < l → l < 1 → 2 < u →
      ∃ G : ℝ → ℝ, ContDiff ℝ ∞ G ∧ EqOn G F phaseInterval ∧
        ∀ p ≤ P, ∀ x ∈ Icc l u,
          ‖iteratedDeriv (p + 1) G x - iteratedDeriv p (modelPhase σ) x‖ ≤ C * δ := by
  obtain ⟨C, hC, hext⟩ := exists_approximateModelPhase_extension P
  refine ⟨C, hC, ?_⟩
  intro δ hδ σ F hF l u hl0 hl hu
  have hlu : l < u := by linarith
  obtain ⟨g, hg, heq, hb⟩ := hext δ hδ σ F hF (l / 2) (2 * u)
    (by positivity) (by linarith) (by linarith)
  obtain ⟨G, hG, hGg⟩ := SmoothExtension.exists_contDiff_eqOn_Icc l u hl0 hlu g hg
  have hJ : UniqueDiffOn ℝ (Icc l u) := uniqueDiffOn_Icc hlu
  have hK : UniqueDiffOn ℝ (Icc (l / 2) (2 * u)) := uniqueDiffOn_Icc (by linarith)
  have hJK : Icc l u ⊆ Icc (l / 2) (2 * u) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hmodel : ContDiffOn ℝ ∞ (modelPhase σ) (Icc (l / 2) (2 * u)) := by
    intro x hx
    exact (Real.contDiffAt_rpow_const_of_ne (ne_of_gt (show 0 < x by linarith
      [hx.1]))).contDiffWithinAt
  refine ⟨G, hG, ?_, ?_⟩
  · intro x hx
    exact (hGg ⟨hl.le.trans hx.1, hx.2.trans hu.le⟩).trans (heq hx)
  · intro p hp x hx
    have hphase := iteratedDerivWithin_eq_of_restriction hJK hJ hK hg hGg.symm hx (p + 1)
    have hreference := iteratedDerivWithin_eq_of_restriction hJK hJ hK hmodel
      (fun _ _ => rfl) hx p
    have hGa : ContDiffAt ℝ (p + 1) G x := hG.contDiffAt.of_le
      (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1))
    have hma : ContDiffAt ℝ p (modelPhase σ) x :=
      Real.contDiffAt_rpow_const_of_ne (ne_of_gt (hl0.trans_le hx.1))
    rw [iteratedDerivWithin_eq_iteratedDeriv hJ hGa hx] at hphase
    rw [iteratedDerivWithin_eq_iteratedDeriv hJ hma hx] at hreference
    rw [← hphase, ← hreference]
    exact hb p hp x (hJK hx)

end
end Expdb
