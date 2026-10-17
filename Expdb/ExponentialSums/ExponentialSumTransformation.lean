module

import Expdb.Basic.AutomaticUniformity
public import Expdb.Basic.DerivativeAsymptotics
public import Expdb.ExponentialSums.DualPhase
public import Expdb.ExponentialSums.SmoothCompletion
public import Expdb.Mathlib.Poisson
public import Expdb.Mathlib.NonstationaryPhase
public import Expdb.Mathlib.QuadraticNormalization
public import Expdb.Mathlib.StationaryPhase

/-!
# The quantitative van der Corput transformation

Proves a quantitative van der Corput transformation for smoothly weighted exponential sums,
combining Poisson summation with stationary and nonstationary phase estimates.
-/

@[expose] public section

namespace Expdb.BProcess

open StationaryPhase

open Filter Topology Set MeasureTheory Complex
open scoped ContDiff FourierTransform SchwartzMap

noncomputable section

section AmplitudeBounds

/-! ## Derivative bounds and convergence along sequences -/

private theorem iteratedDeriv_sub_const_eq (a : ℝ) (k : ℕ) (x : ℝ) :
    iteratedDeriv k (fun y ↦ y - a) x =
      if k = 0 then x - a else if k = 1 then 1 else 0 := by
  by_cases hk : k = 0
  · simp [hk]
  · have h := iteratedDeriv_const_add (f := id) (x := x) (Nat.pos_of_ne_zero hk) (-a)
    simpa [sub_eq_add_neg, add_comm, iteratedDeriv_id, hk] using h

private theorem derivativesClose_sub_const {v z : ℕ → ℝ}
    (hb : VariableObject.IsBounded (fun i ↦ z i - v i)) (n : ℕ) :
    DerivativesClose (fun i y ↦ y - v i) (fun i y ↦ y - v i) z z n := by
  intro k _
  simp only [iteratedDeriv_sub_const_eq]
  rcases k with _ | _ | k
  · exact ⟨IsEqUpToInfinitesimal.refl _, hb⟩
  · exact ⟨IsEqUpToInfinitesimal.refl _, VariableObject.IsBounded.const 1⟩
  · exact ⟨IsEqUpToInfinitesimal.refl _, VariableObject.IsBounded.const 0⟩

section Sequences

variable {H : ℕ → ℝ → ℝ} {c C : ℝ}

/-- Sequence-level bound on the derivatives of the amplitude at points whose preimages stay
bounded. -/
private theorem amplitude_seq_bound (hH : ∀ i, ContDiff ℝ ∞ (H i)) (hc : 0 < c)
    (hcurv : ∀ i x, -C ≤ iteratedDeriv 2 (H i) x ∧ iteratedDeriv 2 (H i) x ≤ -c)
    (hbdd : ∀ k, 2 ≤ k → ∃ K, ∀ᶠ i in atTop, ∀ x, ‖iteratedDeriv k (H i) x‖ ≤ K)
    {w : ℝ → ℝ} (hw : ContDiff ℝ ∞ w) (hwb : ∀ k, ∃ K, ∀ x, ‖iteratedDeriv k w x‖ ≤ K)
    {v y : ℕ → ℝ} {M : ℝ} (hv : ∀ i, |v i| ≤ M)
    (hz : ∀ i, |quadInverse (H i) (v i) (y i)| ≤ M) (k : ℕ) :
    ∃ B, ∀ᶠ i in atTop, ‖iteratedDeriv k
      (fun t ↦ w (quadInverse (H i) (v i) t) * deriv (quadInverse (H i) (v i)) t) (y i)‖ ≤ B := by
  have hC : 0 < C := by linarith [(hcurv 0 0).1, (hcurv 0 0).2]
  set X : ℕ → ℝ → ℝ := fun i ↦ quadInverse (H i) (v i) with hXdef
  set Y : ℕ → ℝ → ℝ := fun i ↦ phaseQuadraticCoordinate (H i) (v i) with hYdef
  have hspec := fun i ↦ quadInverse_spec (hH i) hc (hcurv i) (v i)
  have hX : ∀ i, ContDiff ℝ ∞ (X i) := fun i ↦ (hspec i).1
  have hYX : ∀ i t, Y i (X i t) = t := fun i ↦ (hspec i).2.1
  have hY : ∀ i, ContDiff ℝ ∞ (Y i) := fun i ↦
    contDiff_phaseQuadraticCoordinate (hH i) hc (hcurv i) (v i)
  set z : ℕ → ℝ := fun i ↦ X i (y i) with hzdef
  -- derivative bounds and convergence of the quadratic factor
  set Q : ℕ → ℝ → ℝ := fun i x ↦ -phaseQuadraticFactor (H i) (v i) x with hQdef
  have hQs : ∀ i, ContDiff ℝ ∞ (Q i) := fun i ↦
    (contDiff_phaseQuadraticFactor (H i) (v i) (hH i)).neg
  have hQb : ∀ j, ∃ K, ∀ᶠ i in atTop, ∀ x, ‖iteratedDeriv j (Q i) x‖ ≤ K := by
    intro j
    obtain ⟨K, hK⟩ := hbdd (j + 2) (by omega)
    refine ⟨K, hK.mono fun i hi x ↦ ?_⟩
    simp only [Q, iteratedDeriv_fun_neg, norm_neg]
    exact norm_iteratedDeriv_phaseQuadraticFactor_le (H i) (v i) x j (hH i)
      (fun t _ ↦ hi _)
  have hQlow : ∀ i x, c ≤ Q i x := fun i x ↦ by
    have := (phaseQuadraticFactor_mem (hH i) (hcurv i) (v i) x).2
    simp only [Q]
    linarith
  set S : ℕ → ℝ → ℝ := fun i x ↦ Real.sqrt (Q i x) with hSdef
  have hSs : ∀ i, ContDiff ℝ ∞ (S i) := fun i ↦
    (hQs i).sqrt fun x ↦ (lt_of_lt_of_le hc (hQlow i x)).ne'
  -- boundedness of the preimage points
  have hzb : VariableObject.IsBounded z :=
    ⟨M,
      Eventually.of_forall fun i ↦ by simpa [Real.norm_eq_abs] using hz i⟩
  have hzv : VariableObject.IsBounded (fun i ↦ z i - v i) :=
    ⟨2 * M, Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs]
      have := hz i
      have := hv i
      simp only [z, X] at *
      calc |quadInverse (H i) (v i) (y i) - v i|
          ≤ |quadInverse (H i) (v i) (y i)| + |v i| := abs_sub _ _
        _ ≤ 2 * M := by linarith⟩
  have hYj : ∀ n, DerivativesClose Y Y z z n := by
    intro n
    have hSj : DerivativesClose S S z z n := by
      refine DerivativesClose.of_mul_self_eq hSs hSs (fun _ _ ↦ Real.sqrt_nonneg _)
        (fun _ _ ↦ Real.sqrt_nonneg _)
        (fun i x ↦ Real.mul_self_sqrt (hc.le.trans (hQlow i x)))
        (fun i x ↦ Real.mul_self_sqrt (hc.le.trans (hQlow i x)))
        (DerivativesClose.of_uniform_bound fun j _ ↦ hQb j) (Real.sqrt_pos.2 hc)
        (Eventually.of_forall fun i ↦ Real.sqrt_le_sqrt (hQlow i (z i)))
    exact DerivativesClose.mul (f := fun i y ↦ y - v i) (g := S)
      (f' := fun i y ↦ y - v i) (g' := S)
      (fun _ ↦ contDiff_id.sub contDiff_const) hSs (fun _ ↦ contDiff_id.sub contDiff_const)
      hSs (derivativesClose_sub_const hzv n) hSj
  -- derivative bounds and convergence of the inverse
  have hXj : ∀ n, DerivativesClose X X y y (n + 1) := by
    intro n
    refine DerivativesClose.inverse hX hX hY hY hYX hYX (hYj (n + 1)).shift
      (IsEqUpToInfinitesimal.refl _) hzb
      (div_pos hc (Real.sqrt_pos.2 hC)) (Eventually.of_forall fun i ↦ ?_)
    have := le_deriv_phaseQuadraticCoordinate (hH i) hc (hcurv i) (v i) (z i)
    rw [Real.norm_of_nonneg ((div_pos hc (Real.sqrt_pos.2 hC)).le.trans this)]
    exact this
  -- derivative bounds and convergence of the amplitude
  have hwj : DerivativesClose (fun _ ↦ w) (fun _ ↦ w) z z k := DerivativesClose.of_uniform_bound
    fun j _ ↦ by
      obtain ⟨K, hK⟩ := hwb j
      exact ⟨K, Eventually.of_forall fun _ ↦ hK⟩
  have hcomp := DerivativesClose.comp (g := fun _ ↦ w) (g' := fun _ ↦ w) (f := X) (f' := X)
    (x := y) (x' := y) (fun _ ↦ hw) (fun _ ↦ hw) hX hX hwj
    ((hXj k).mono (by omega)).positiveDerivativesClose
  have hamp := DerivativesClose.mul (f := fun i ↦ w ∘ X i) (g := fun i ↦ deriv (X i))
    (f' := fun i ↦ w ∘ X i) (g' := fun i ↦ deriv (X i))
    (fun i ↦ hw.comp (hX i)) (fun i ↦ (hX i).deriv' (n := ∞))
    (fun i ↦ hw.comp (hX i)) (fun i ↦ (hX i).deriv' (n := ∞)) hcomp (hXj k).shift
  exact (hamp k le_rfl).2

/-- Uniform bounds on the derivatives of the stationary-phase amplitude. -/
theorem exists_amplitude_bound (hH : ∀ i, ContDiff ℝ ∞ (H i)) (hc : 0 < c)
    (hcurv : ∀ i x, -C ≤ iteratedDeriv 2 (H i) x ∧ iteratedDeriv 2 (H i) x ≤ -c)
    (hbdd : ∀ k, 2 ≤ k → ∃ K, ∀ᶠ i in atTop, ∀ x, ‖iteratedDeriv k (H i) x‖ ≤ K)
    {w : ℝ → ℝ} (hw : ContDiff ℝ ∞ w) {p q : ℝ} (hpq : p ≤ q)
    (hw0 : ∀ x, x ∉ Icc p q → w x = 0) {M : ℝ} (k : ℕ) :
    ∃ B, ∀ᶠ i in atTop, ∀ v, |v| ≤ M → ∀ t, ‖iteratedDeriv k
      (fun t ↦ w (quadInverse (H i) v t) * deriv (quadInverse (H i) v) t) t‖ ≤ B := by
  classical
  have hwc : HasCompactSupport w := by
    apply HasCompactSupport.intro (isCompact_Icc (a := p) (b := q))
    exact hw0
  have hwb : ∀ k, ∃ K, ∀ x, ‖iteratedDeriv k w x‖ ≤ K := fun k ↦
    exists_bound_iteratedDeriv_of_hasCompactSupport hw hwc k
  by_cases hM : 0 ≤ M
  swap
  · exact ⟨0, Eventually.of_forall fun i v hv ↦ absurd ((abs_nonneg v).trans hv) hM⟩
  set M' : ℝ := max M (max |p - 1| |q + 1|) with hM'
  let P := {v : ℝ // |v| ≤ M} × ℝ
  haveI : Nonempty P := ⟨(⟨0, by simpa using hM⟩, 0)⟩
  let g : ℕ → P → ℝ := fun i r ↦ iteratedDeriv k
    (fun t ↦ w (quadInverse (H i) r.1 t) * deriv (quadInverse (H i) r.1) t) r.2
  suffices hsuff : ∃ B, ∀ᶠ i in atTop, ∀ r, ‖g i r‖ ≤ B by
    obtain ⟨B, hB⟩ := hsuff
    exact ⟨B, hB.mono fun i hi v hv t ↦ hi (⟨v, hv⟩, t)⟩
  let E : VariableObject (Set P) := fun _ ↦ Set.univ
  let f : VariableFunction (fun i ↦ E i) ℝ := fun i q ↦ g i q
  suffices hf : f.IsChoicewiseBounded by
    obtain ⟨C, hC⟩ :=
      (f.isChoicewiseBounded_iff_eventually_uniform E (fun _ ↦ Set.univ_nonempty)).1 hf
    exact ⟨C, hC.mono fun i hi q ↦ hi ⟨q, Set.mem_univ q⟩⟩
  intro choice
  let r : ℕ → P := fun i ↦ (choice i).val
  change ∃ B, ∀ᶠ i in atTop, ‖g i (r i)‖ ≤ B
  set v : ℕ → ℝ := fun i ↦ (r i).1.1
  set y : ℕ → ℝ := fun i ↦ (r i).2
  -- replace bad points
  set good : ℕ → Prop := fun i ↦ quadInverse (H i) (v i) (y i) ∈ Icc (p - 1) (q + 1)
  set y' : ℕ → ℝ := fun i ↦ if good i then y i else phaseQuadraticCoordinate (H i) (v i) p
  have hspec := fun i ↦ quadInverse_spec (hH i) hc (hcurv i) (v i)
  have hz' : ∀ i, |quadInverse (H i) (v i) (y' i)| ≤ M' := by
    intro i
    by_cases hg : good i
    · simp only [y', if_pos hg]
      have := hg
      simp only [good, mem_Icc] at this
      rw [abs_le]
      constructor
      · have : -|p - 1| ≤ p - 1 := neg_abs_le _
        have : |p - 1| ≤ M' := le_max_of_le_right (le_max_left _ _)
        linarith
      · have : q + 1 ≤ |q + 1| := le_abs_self _
        have : |q + 1| ≤ M' := le_max_of_le_right (le_max_right _ _)
        linarith
    · simp only [y', if_neg hg, (hspec i).2.2]
      rw [abs_le]
      constructor
      · have : -|p - 1| ≤ p - 1 := neg_abs_le _
        have : |p - 1| ≤ M' := le_max_of_le_right (le_max_left _ _)
        linarith
      · have : q + 1 ≤ |q + 1| := le_abs_self _
        have : |q + 1| ≤ M' := le_max_of_le_right (le_max_right _ _)
        linarith
  have hv' : ∀ i, |v i| ≤ M' := fun i ↦ (r i).1.2.trans (le_max_left _ _)
  obtain ⟨B, hB⟩ := amplitude_seq_bound hH hc hcurv hbdd hw hwb hv' hz' k
  refine ⟨B, hB.mono fun i hi ↦ ?_⟩
  by_cases hg : good i
  · simpa [g, y', if_pos hg] using hi
  · have hzero : iteratedDeriv k
        (fun t ↦ w (quadInverse (H i) (v i) t) * deriv (quadInverse (H i) (v i)) t) (y i) = 0 := by
      have hcont : Continuous (quadInverse (H i) (v i)) := (hspec i).1.continuous
      have hopen : IsOpen {t | quadInverse (H i) (v i) t ∉ Icc p q} :=
        (isClosed_Icc.preimage hcont).isOpen_compl
      have hmem : y i ∈ {t | quadInverse (H i) (v i) t ∉ Icc p q} := by
        intro hmem
        apply hg
        exact ⟨by linarith [hmem.1], by linarith [hmem.2]⟩
      have hev : (fun t ↦ w (quadInverse (H i) (v i) t) * deriv (quadInverse (H i) (v i)) t)
          =ᶠ[𝓝 (y i)] fun _ ↦ 0 := by
        filter_upwards [hopen.mem_nhds hmem] with t ht
        simp [hw0 _ ht]
      rw [hev.iteratedDeriv_eq]
      simp
    have : g i (r i) = 0 := hzero
    rw [this, norm_zero]
    exact (norm_nonneg _).trans hi

end Sequences

end AmplitudeBounds

section SmoothPoisson

/-- Poisson summation for a smoothly weighted exponential sum. -/
private theorem smooth_sum_eq_poisson {φ : ℝ → ℝ} (hφ : ContDiff ℝ ∞ φ) {T N : ℝ} (hT : 0 < T)
    (hN : 0 < N) {B : ℕ} (hB : 9 / 4 * N ≤ B) :
    ∑ n ∈ Finset.Icc 0 B, (phaseCutoff ((n : ℝ) / N) : ℂ) * 𝐞 (T * φ ((n : ℝ) / N)) =
      N * ∑' m : ℤ, ∫ u : ℝ, (phaseCutoff u : ℂ) * 𝐞 (T * (φ u - N * m / T * u)) := by
  set g : ℝ → ℂ := fun u ↦ (phaseCutoff u : ℂ) * 𝐞 (T * φ u) with hg
  have hchar : ContDiff ℝ ∞ (fun u ↦ (𝐞 (T * φ u) : ℂ)) := by
    simp only [Real.fourierChar_apply]
    exact Complex.contDiff_exp.comp
      ((Complex.ofRealCLM.contDiff.comp (contDiff_const.mul
        (contDiff_const.mul hφ))).mul contDiff_const)
  have hgs : ContDiff ℝ ∞ g :=
    (Complex.ofRealCLM.contDiff.comp phaseCutoff.contDiff).mul hchar
  have hgc : HasCompactSupport g := by
    apply HasCompactSupport.mul_right
    exact phaseCutoff.hasCompactSupport.comp_left Complex.ofReal_zero
  let f : 𝓢(ℝ, ℂ) := hgc.toSchwartzMap hgs
  have hf : ∀ u, f u = g u := fun u ↦ by simp [f]
  have hP := poisson_scaled f hN
  simp only [hf] at hP
  -- the left side
  have hleft : (∑' n : ℤ, g ((n : ℝ) / N)) =
      ∑ n ∈ Finset.Icc 0 B, (phaseCutoff ((n : ℝ) / N) : ℂ) * 𝐞 (T * φ ((n : ℝ) / N)) := by
    rw [tsum_eq_sum (s := (Finset.Icc 0 B).image (fun n : ℕ ↦ (n : ℤ)))]
    · rw [Finset.sum_image (fun a _ b _ h ↦ by exact_mod_cast h)]
      rfl
    · intro n hn
      have hz : phaseCutoff ((n : ℝ) / N) = 0 := by
        apply phaseCutoff_eq_zero
        rintro ⟨h1, h2⟩
        apply hn
        rw [Finset.mem_image]
        have hn0 : 0 ≤ n := by
          by_contra h
          push Not at h
          have : (n : ℝ) / N < 0 := div_neg_of_neg_of_pos (by exact_mod_cast h) hN
          linarith
        refine ⟨n.toNat, Finset.mem_Icc.2 ⟨Nat.zero_le _, ?_⟩, Int.toNat_of_nonneg hn0⟩
        have : (n : ℝ) < 9 / 4 * N := by rwa [div_lt_iff₀ hN] at h2
        have : (n : ℝ) < B := lt_of_lt_of_le this hB
        have : n < (B : ℤ) := by exact_mod_cast this
        omega
      simp [g, hz]
  rw [← hleft, hP, Complex.real_smul]
  congr 1
  apply tsum_congr
  intro m
  apply integral_congr_ae
  filter_upwards with u
  simp only [g, Circle.smul_def, smul_eq_mul]
  rw [mul_left_comm, ← Circle.coe_mul, ← AddChar.map_add_eq_mul]
  congr 3
  field_simp
  ring

end SmoothPoisson

section Transformation

/-- The amplitude of the stationary main term at frequency `ν`. -/
def mainAmplitude (H : ℝ → ℝ) (ν : ℝ) : ℝ :=
  phaseCutoff (phaseInverse H ν) / Real.sqrt (-iteratedDeriv 2 H (phaseInverse H ν))

/-- The stationary main term at frequency `ν`. -/
def mainTerm (H : ℝ → ℝ) (T ν : ℝ) : ℂ :=
  𝐞 (T * legendre H ν) * (Complex.exp (-(Real.pi : ℂ) / 4 * I) *
    (mainAmplitude H ν : ℂ) / (Real.sqrt T : ℂ))

/-- The transformed stationary-phase amplitude centred at `v`. -/
def stationaryAmplitude (H : ℝ → ℝ) (v : ℝ) (y : ℝ) : ℝ :=
  phaseCutoff (quadInverse H v y) * deriv (quadInverse H v) y

private theorem summable_inv_one_add_abs_sq : Summable (fun m : ℤ ↦ 1 / (1 + |(m : ℝ)|) ^ 2) := by
  have h := (Real.summable_one_div_nat_add_rpow 1 2).2 (by norm_num)
  rw [summable_int_iff_summable_nat_and_neg]
  constructor <;> simpa [Int.cast_natCast, Int.cast_neg, abs_neg, Real.rpow_two,
    abs_of_nonneg, add_comm] using h

/-- A common bound for the finitely many cutoff derivatives used in the frequency estimates. -/
theorem exists_phaseCutoff_bound (n : ℕ) :
    ∃ W, ∀ k ≤ n, ∀ x, |iteratedDeriv k (phaseCutoff : ℝ → ℝ) x| ≤ W := by
  classical
  choose B hB using fun k : ℕ ↦ exists_bound_iteratedDeriv_of_hasCompactSupport
    phaseCutoff.contDiff phaseCutoff.hasCompactSupport k
  refine ⟨∑ k ∈ Finset.range (n + 1), B k, fun k hk x ↦ ?_⟩
  calc
    |iteratedDeriv k (phaseCutoff : ℝ → ℝ) x| ≤ B k := by
      simpa only [Real.norm_eq_abs] using hB k x
    _ ≤ _ := Finset.single_le_sum (fun j _ ↦ (norm_nonneg _).trans (hB j 0))
      (Finset.mem_range.2 (by omega))


section

variable {H : ℝ → ℝ} {c C : ℝ}

/-- The stationary estimate at a single frequency. Bounds on the second and fourth
derivatives of the transformed amplitude give an error of order `T⁻³ᐟ²`. -/
private theorem norm_integral_sub_mainTerm_le (hH : ContDiff ℝ ∞ H) (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    {T ν M0 Bamp : ℝ} (hT : 0 < T) (hv : |phaseInverse H ν| ≤ M0)
    (hamp : ∀ v, |v| ≤ M0 → ∀ y, ‖iteratedDeriv 2 (stationaryAmplitude H v) y‖ ≤ Bamp ∧
      ‖iteratedDeriv 4 (stationaryAmplitude H v) y‖ ≤ Bamp) :
    ‖(∫ u : ℝ, (phaseCutoff u : ℂ) * 𝐞 (T * (H u - ν * u))) - mainTerm H T ν‖ ≤
      4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp / (T * Real.sqrt T) := by
  set v := phaseInverse H ν with hvdef
  obtain ⟨-, -, hfg⟩ := phaseInverse_spec hH hc (fun x ↦ (hcurv x).2)
  have hνv : deriv H v = ν := hfg ν
  obtain ⟨hX, hYX, -⟩ := quadInverse_spec hH hc hcurv v
  have hM0 : 0 ≤ M0 := (abs_nonneg _).trans hv
  have hC : 0 ≤ C := by linarith [(hcurv 0).1, (hcurv 0).2]
  have hsupp : ∀ y, (9 / 4 + M0) * Real.sqrt C < |y| → phaseCutoff (quadInverse H v y) = 0 := by
    intro y hy
    by_contra hne
    have hmem : quadInverse H v y ∈ Ioo (3 / 4 : ℝ) (9 / 4) := by
      by_contra h
      exact hne (phaseCutoff_eq_zero h)
    have hb := abs_phaseQuadraticCoordinate_le hH hcurv v (quadInverse H v y)
    rw [hYX] at hb
    have h1 : |quadInverse H v y - v| ≤ 9 / 4 + M0 := by
      calc |quadInverse H v y - v| ≤ |quadInverse H v y| + |v| := abs_sub _ _
        _ ≤ 9 / 4 + M0 := by
          rw [abs_of_pos (by linarith [hmem.1])]
          linarith [hmem.2]
    have := mul_le_mul_of_nonneg_right h1 (Real.sqrt_nonneg C)
    linarith
  have hst := StationaryPhase.norm_integral_stationary_sub_le
    (w := phaseCutoff) (X := quadInverse H v)
    (lam := T) (R := (9 / 4 + M0) * Real.sqrt C) (B := Bamp) hH hc hcurv phaseCutoff.contDiff
    hX hYX hT (by positivity) hsupp (fun y ↦ (hamp v hv y).1) (fun y ↦ (hamp v hv y).2)
  rw [hνv] at hst
  convert hst using 3
  simp only [mainTerm, mainAmplitude, legendre, ← hvdef]

/-- The nonstationary estimate at a single frequency. The cutoff is supported in
`(3/4, 9/4)`; the larger interval `[5/8, 19/8]` allows integration by parts without
boundary terms. Bounds on the inverse phase derivative and its first two derivatives
give a summable bound via `Expdb.Mathlib.NonstationaryPhase`. -/
private theorem norm_integral_nonstationary_freq_le (hH : ContDiff ℝ ∞ H) (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    {T ν δ C3 W : ℝ} (hT : 0 < T) (hδ : 0 < δ)
    (hsep : ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8), δ * (1 + |ν|) ≤ |deriv H x - ν|)
    (h3 : ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8), |iteratedDeriv 3 H x| ≤ C3)
    (hW : ∀ k ≤ 2, ∀ x, |iteratedDeriv k (phaseCutoff : ℝ → ℝ) x| ≤ W) :
    ‖∫ u : ℝ, (phaseCutoff u : ℂ) * 𝐞 (T * (H u - ν * u))‖ ≤
      (7 / 4 * W * (1 + 3 * C / δ + (C3 / δ + 3 * C ^ 2 / δ ^ 2))) /
        (T ^ 2 * (δ * (1 + |ν|)) ^ 2) := by
  have hd : 0 < δ * (1 + |ν|) := by positivity
  have hdδ : δ ≤ δ * (1 + |ν|) := by nlinarith [abs_nonneg ν]
  have hC2 : ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8), |iteratedDeriv 2 H x| ≤ C := fun x _ ↦ by
    rw [abs_le]; constructor <;> linarith [(hcurv x).1, (hcurv x).2]
  have hC : 0 ≤ C := (abs_nonneg _).trans (hC2 1 ⟨by norm_num, by norm_num⟩)
  have hC3 : 0 ≤ C3 := (abs_nonneg _).trans (h3 1 ⟨by norm_num, by norm_num⟩)
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0 (by norm_num) 0)
  have h := NonstationaryPhase.norm_integral_nonstationary_le
    (a := 5 / 8) (a' := 3 / 4) (b' := 9 / 4) (b := 19 / 8)
    (by norm_num) (by norm_num) (by norm_num) hH phaseCutoff.contDiff
    (fun x hx ↦ phaseCutoff_eq_zero fun h ↦ hx ⟨h.1.le, h.2.le⟩) hT hd hsep hC2 h3 hW
  refine h.trans ?_
  apply div_le_div_of_nonneg_right _ (by positivity)
  have e1 : 3 * C / (δ * (1 + |ν|)) ≤ 3 * C / δ :=
    div_le_div_of_nonneg_left (by positivity) hδ hdδ
  have e2 : C3 / (δ * (1 + |ν|)) ≤ C3 / δ := div_le_div_of_nonneg_left hC3 hδ hdδ
  have e3 : 3 * C ^ 2 / (δ * (1 + |ν|)) ^ 2 ≤ 3 * C ^ 2 / δ ^ 2 :=
    div_le_div_of_nonneg_left (by positivity) (by positivity)
      (pow_le_pow_left₀ hδ.le hdδ 2)
  norm_num
  gcongr


/-- Poisson summation followed by stationary phase at the frequencies in `[P₀, P₁)` and
nonstationary phase at all the others. -/
theorem norm_smooth_sum_sub_main_le (hH : ContDiff ℝ ∞ H) (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    {T N : ℝ} (hT : 0 < T) (hN : 0 < N) (hNT : N ≤ T) {B : ℕ} (hB : 9 / 4 * N ≤ B)
    {P₀ P₁ : ℕ} {M0 Bamp δ C3 W : ℝ} (hδ : 0 < δ)
    (hstat : ∀ m ∈ Finset.Ico P₀ P₁, |phaseInverse H (N * m / T)| ≤ M0)
    (hamp : ∀ v, |v| ≤ M0 → ∀ y, ‖iteratedDeriv 2 (stationaryAmplitude H v) y‖ ≤ Bamp ∧
      ‖iteratedDeriv 4 (stationaryAmplitude H v) y‖ ≤ Bamp)
    (hsep : ∀ m : ℤ, (m < P₀ ∨ (P₁ : ℤ) ≤ m) → ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8),
      δ * (1 + |N * m / T|) ≤ |deriv H x - N * m / T|)
    (h3 : ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8), |iteratedDeriv 3 H x| ≤ C3)
    (hW : ∀ k ≤ 2, ∀ x, |iteratedDeriv k (phaseCutoff : ℝ → ℝ) x| ≤ W) :
    ‖∑ n ∈ Finset.Icc 0 B, (phaseCutoff ((n : ℝ) / N) : ℂ) * 𝐞 (T * H ((n : ℝ) / N)) -
        N * ∑ m ∈ Finset.Ico P₀ P₁, mainTerm H T (N * m / T)‖ ≤
      N * (((P₁ - P₀ : ℕ) : ℝ) *
          (4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp / (T * Real.sqrt T))) +
        (7 / 4 * W * (1 + 3 * C / δ + (C3 / δ + 3 * C ^ 2 / δ ^ 2))) *
          (∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2) / (δ ^ 2 * N) := by
  classical
  set Es : ℝ := 4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp / (T * Real.sqrt T)
    with hEs
  set Kn : ℝ := 7 / 4 * W * (1 + 3 * C / δ + (C3 / δ + 3 * C ^ 2 / δ ^ 2)) with hKn
  set S0 : ℝ := ∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2 with hS0
  rw [smooth_sum_eq_poisson hH hT hN hB]
  set I : ℤ → ℂ := fun m ↦ ∫ u : ℝ, (phaseCutoff u : ℂ) * 𝐞 (T * (H u - N * m / T * u))
    with hI
  set S : Finset ℤ := (Finset.Ico P₀ P₁).image (fun n : ℕ ↦ (n : ℤ)) with hS
  set f : ℤ → ℂ := fun m ↦ if m ∈ S then mainTerm H T (N * m / T) else 0 with hf
  set E : ℤ → ℂ := fun m ↦ I m - f m with hE
  set g : ℤ → ℝ := fun m ↦ (if m ∈ S then Es else 0) +
    Kn / (δ ^ 2 * N ^ 2) * (1 / (1 + |(m : ℝ)|) ^ 2) with hg
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0 (by norm_num) 0)
  have hC : 0 ≤ C := by linarith [(hcurv 0).1, (hcurv 0).2]
  have hC3 : 0 ≤ C3 := (abs_nonneg _).trans (h3 1 ⟨by norm_num, by norm_num⟩)
  have hKn : 0 ≤ Kn := by positivity
  -- pointwise bounds
  have hmemS : ∀ m : ℤ, m ∈ S ↔ ∃ n ∈ Finset.Ico P₀ P₁, (n : ℤ) = m := by
    intro m; simp [S]
  have hbound : ∀ m, ‖E m‖ ≤ g m := by
    intro m
    have h2 : 0 ≤ Kn / (δ ^ 2 * N ^ 2) * (1 / (1 + |(m : ℝ)|) ^ 2) := by positivity
    by_cases hm : m ∈ S
    · obtain ⟨n, hn, rfl⟩ := (hmemS m).1 hm
      simp only [E, f, g, if_pos hm]
      have := norm_integral_sub_mainTerm_le hH hc hcurv (ν := N * (n : ℤ) / T) hT
        (by simpa using hstat n hn) hamp
      simp only [I]
      linarith
    · simp only [E, f, g, if_neg hm, sub_zero, zero_add]
      have hout : m < P₀ ∨ (P₁ : ℤ) ≤ m := by
        by_contra hcon
        push Not at hcon
        apply hm
        rw [hmemS]
        refine ⟨m.toNat, Finset.mem_Ico.2 ⟨by omega, by omega⟩, by omega⟩
      have hns := norm_integral_nonstationary_freq_le hH hc hcurv (ν := N * m / T) hT hδ
        (hsep m hout) h3 hW
      refine hns.trans ?_
      have hden : δ ^ 2 * N ^ 2 * (1 + |(m : ℝ)|) ^ 2 ≤ T ^ 2 * (δ * (1 + |N * m / T|)) ^ 2 := by
        have h1 : N * (1 + |(m : ℝ)|) ≤ T * (1 + |N * m / T|) := by
          rw [abs_div, abs_mul, abs_of_pos hN, abs_of_pos hT]
          field_simp
          nlinarith [abs_nonneg (m : ℝ)]
        have h2 : (N * (1 + |(m : ℝ)|)) ^ 2 ≤ (T * (1 + |N * m / T|)) ^ 2 :=
          pow_le_pow_left₀ (by positivity) h1 2
        nlinarith [sq_nonneg δ]
      calc Kn / (T ^ 2 * (δ * (1 + |N * m / T|)) ^ 2)
          ≤ Kn / (δ ^ 2 * N ^ 2 * (1 + |(m : ℝ)|) ^ 2) :=
            div_le_div_of_nonneg_left hKn (by positivity) hden
        _ = Kn / (δ ^ 2 * N ^ 2) * (1 / (1 + |(m : ℝ)|) ^ 2) := by
            field_simp
  have hgs : Summable g := by
    refine Summable.add ?_ (summable_inv_one_add_abs_sq.mul_left _)
    exact summable_of_ne_finset_zero (s := S) (fun m hm ↦ by simp [hm])
  have hEs' : Summable E := Summable.of_norm_bounded hgs hbound
  have hfs : Summable f := summable_of_ne_finset_zero (s := S) (fun m hm ↦ by simp [f, hm])
  have hIf : ∀ m, I m = f m + E m := fun m ↦ by simp [E]
  have hsumf : ∑' m, f m = ∑ m ∈ Finset.Ico P₀ P₁, mainTerm H T (N * m / T) := by
    rw [tsum_eq_sum (s := S) (fun m hm ↦ by simp [f, hm])]
    simp only [f, Finset.sum_ite_mem, Finset.inter_self]
    rw [hS, Finset.sum_image (fun a _ b _ h ↦ by exact_mod_cast h)]
    simp only [Int.cast_natCast]
  have htsum : ∑' m, I m = ∑ m ∈ Finset.Ico P₀ P₁, mainTerm H T (N * m / T) + ∑' m, E m := by
    rw [tsum_congr hIf, hfs.tsum_add hEs', hsumf]
  have hcard : S.card = P₁ - P₀ := by
    rw [hS, Finset.card_image_of_injective _ (fun a b h ↦ by exact_mod_cast h), Nat.card_Ico]
  have hsumg : ∑' m, g m = ((P₁ - P₀ : ℕ) : ℝ) * Es + Kn / (δ ^ 2 * N ^ 2) * S0 := by
    rw [(summable_of_ne_finset_zero (s := S) (fun m hm ↦ by simp [hm])).tsum_add
      (summable_inv_one_add_abs_sq.mul_left _), tsum_mul_left,
      tsum_eq_sum (s := S) (fun m hm ↦ by simp [hm])]
    simp only [Finset.sum_ite_mem, Finset.inter_self, Finset.sum_const, hcard, nsmul_eq_mul, S0]
  have hnormE : ‖∑' m, E m‖ ≤ ((P₁ - P₀ : ℕ) : ℝ) * Es + Kn / (δ ^ 2 * N ^ 2) * S0 := by
    rw [← hsumg]
    exact tsum_of_norm_bounded hgs.hasSum hbound
  rw [htsum, mul_add, add_sub_cancel_left, norm_mul, Complex.norm_real, Real.norm_of_nonneg hN.le]
  calc N * ‖∑' m, E m‖ ≤ N * (((P₁ - P₀ : ℕ) : ℝ) * Es + Kn / (δ ^ 2 * N ^ 2) * S0) :=
        mul_le_mul_of_nonneg_left hnormE hN.le
    _ = _ := by field_simp

end

end Transformation

end

end Expdb.BProcess
