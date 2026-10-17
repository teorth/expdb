module

import Mathlib.Data.Nat.Log

public import Mathlib.Analysis.Complex.Circle
public import Mathlib.Order.Interval.Finset.Nat

import Mathlib.Algebra.Order.Chebyshev

/-!
# Van der Corput differencing

Adapted from `MoltResearch/Discrepancy/ExpSums.lean` in
[ProofFleet/prooffleet](https://github.com/ProofFleet/prooffleet) at commit
`9f6f78688f83561d04616ab88e93b9cf81531e47`, licensed under Apache-2.0
(see the repository's `LICENSE`).

Adaptations for Expdb: retained differencing and its dependencies, made implementation
helpers private, and expressed character identities using Mathlib's Fourier character.
The public interface uses closed natural-number intervals; namespaces and proofs
were updated for the current Mathlib API.
-/

@[expose] public section

namespace Expdb.VanDerCorput

open Complex
open scoped FourierTransform

noncomputable section

private abbrev e (x : ℝ) : ℂ := 𝐞 x

private theorem norm_e (x : ℝ) : ‖e x‖ = 1 := by simp [e]

private theorem e_add (x y : ℝ) : e (x + y) = e x * e y := by
  simp [e, AddChar.map_add_eq_mul]

private theorem e_conj (x : ℝ) : (starRingEnd ℂ) (e x) = e (-x) := by simp [e]

/-- Products against conjugates subtract phases. -/
private theorem e_mul_conj (x y : ℝ) :
    e x * (starRingEnd ℂ) (e y) = e (x - y) := by
  rw [e_conj, ← e_add]
  congr 1

/-- The shift identity: summing a shifted phase over a shifted window. -/
private theorem sum_shift (φ : ℕ → ℝ) {M N : ℕ} (h : ℕ) (hhM : h ≤ M)
    (hMN : M ≤ N + 1) :
    ∑ n ∈ Finset.Ico M (N + 1), e (φ n)
      = ∑ n ∈ Finset.Ico (M - h) (N + 1 - h), e (φ (n + h)) := by
  have hmap : Finset.Ico M (N + 1)
      = (Finset.Ico (M - h) (N + 1 - h)).map
          (addRightEmbedding h) := by
    rw [Finset.map_add_right_Ico]
    congr 1 <;> omega
  rw [hmap, Finset.sum_map]
  simp only [addRightEmbedding_apply]

/-- Discrete Cauchy–Schwarz for complex sums. -/
private theorem norm_sum_sq_le_card_mul (W : Finset ℕ) (f : ℕ → ℂ) :
    ‖∑ n ∈ W, f n‖ ^ 2 ≤ (W.card : ℝ) * ∑ n ∈ W, ‖f n‖ ^ 2 := by
  have h1 : ‖∑ n ∈ W, f n‖ ≤ ∑ n ∈ W, ‖f n‖ := norm_sum_le _ _
  have h2 : (∑ n ∈ W, ‖f n‖) ^ 2 ≤ (W.card : ℝ) * ∑ n ∈ W, ‖f n‖ ^ 2 := by
    exact_mod_cast sq_sum_le_card_mul_sum_sq (s := W) (f := fun n => ‖f n‖)
  calc ‖∑ n ∈ W, f n‖ ^ 2
      ≤ (∑ n ∈ W, ‖f n‖) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ (W.card : ℝ) * ∑ n ∈ W, ‖f n‖ ^ 2 := h2

/-- Expanding the square of an inner sum into conjugate pairs. -/
private theorem sum_norm_sq_expand (W K : Finset ℕ) (c : ℕ → ℕ → ℂ) :
    ∑ n ∈ W, ‖∑ h ∈ K, c n h‖ ^ 2
      ≤ ∑ p ∈ K ×ˢ K, ‖∑ n ∈ W, c n p.1 * (starRingEnd ℂ) (c n p.2)‖ := by
  have hsq : ∀ z : ℂ, (‖z‖ : ℝ) ^ 2 = (z * (starRingEnd ℂ) z).re := by
    intro z
    rw [Complex.mul_conj]
    simp [Complex.sq_norm]
  have hexp : ∀ n, (∑ h ∈ K, c n h) * (starRingEnd ℂ) (∑ h ∈ K, c n h)
      = ∑ p ∈ K ×ˢ K, c n p.1 * (starRingEnd ℂ) (c n p.2) := by
    intro n
    rw [map_sum, Finset.sum_mul_sum, ← Finset.sum_product']
  calc ∑ n ∈ W, ‖∑ h ∈ K, c n h‖ ^ 2
      = ∑ n ∈ W, ((∑ h ∈ K, c n h)
          * (starRingEnd ℂ) (∑ h ∈ K, c n h)).re := by
        refine Finset.sum_congr rfl fun n _ => hsq _
    _ = (∑ n ∈ W, ∑ p ∈ K ×ˢ K,
          c n p.1 * (starRingEnd ℂ) (c n p.2)).re := by
        rw [Complex.re_sum]
        refine Finset.sum_congr rfl fun n _ => ?_
        rw [hexp]
    _ = (∑ p ∈ K ×ˢ K, ∑ n ∈ W,
          c n p.1 * (starRingEnd ℂ) (c n p.2)).re := by
        rw [Finset.sum_comm]
    _ ≤ ‖∑ p ∈ K ×ˢ K, ∑ n ∈ W,
          c n p.1 * (starRingEnd ℂ) (c n p.2)‖ := Complex.re_le_norm _
    _ ≤ ∑ p ∈ K ×ˢ K, ‖∑ n ∈ W,
          c n p.1 * (starRingEnd ℂ) (c n p.2)‖ := norm_sum_le _ _

/-- The off-diagonal fibering: pairs at difference `g` number at most `H`. -/
private theorem sum_offdiag_le {H : ℕ} (F : ℕ → ℝ) (hF : ∀ g, 0 ≤ F g) :
    ∑ p ∈ ((Finset.range H) ×ˢ (Finset.range H)).filter
        (fun p => p.2 < p.1), F (p.1 - p.2)
      ≤ (H : ℝ) * ∑ g ∈ Finset.Ico 1 H, F g := by
  classical
  have hmaps : ∀ p ∈ ((Finset.range H) ×ˢ (Finset.range H)).filter
      (fun p => p.2 < p.1), p.1 - p.2 ∈ Finset.Ico 1 H := by
    intro p hp
    rw [Finset.mem_filter, Finset.mem_product, Finset.mem_range,
      Finset.mem_range] at hp
    rw [Finset.mem_Ico]
    omega
  have hfib := Finset.sum_fiberwise_of_maps_to
    (g := fun p : ℕ × ℕ => p.1 - p.2) (t := Finset.Ico 1 H) hmaps
    (fun p : ℕ × ℕ => F (p.1 - p.2))
  rw [← hfib]
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun g hg => ?_
  have hcong : ∀ p ∈ (((Finset.range H) ×ˢ (Finset.range H)).filter
      (fun p => p.2 < p.1)).filter (fun p => p.1 - p.2 = g),
      F (p.1 - p.2) = F g := by
    intro p hp
    rw [Finset.mem_filter] at hp
    rw [hp.2]
  rw [Finset.sum_congr rfl hcong, Finset.sum_const, nsmul_eq_mul]
  refine mul_le_mul_of_nonneg_right ?_ (hF g)
  have hcard : ((((Finset.range H) ×ˢ (Finset.range H)).filter
      (fun p => p.2 < p.1)).filter (fun p => p.1 - p.2 = g)).card
      ≤ H := by
    have hinj : Set.MapsTo Prod.snd
        (↑((((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => p.2 < p.1)).filter (fun p => p.1 - p.2 = g))
          : Set (ℕ × ℕ))
        ((Finset.range H : Finset ℕ) : Set ℕ) := by
      intro p hp
      rw [Finset.mem_coe, Finset.mem_filter, Finset.mem_filter,
        Finset.mem_product] at hp
      rw [Finset.mem_coe]
      exact hp.1.1.2
    have h2 := Finset.card_le_card_of_injOn Prod.snd hinj ?_
    · rw [Finset.card_range] at h2
      exact h2
    · intro p hp q hq hpq
      rw [Finset.mem_coe, Finset.mem_filter, Finset.mem_filter] at hp hq
      have h3 := hp.2
      have h4 := hq.2
      have h5 := hp.1.2
      have h6 := hq.1.2
      refine Prod.ext ?_ hpq
      omega
  exact_mod_cast hcard

/-- The trivial character value. -/
private theorem e_zero : e 0 = 1 := by
  rw [e]
  norm_num

/-- The off-diagonal pair sum collapses to a difference-phase sum. -/
private theorem weyl_offdiag {φ : ℕ → ℝ} {M N H : ℕ} (hHM : H ≤ M) (hMN : M ≤ N)
    {h h' : ℕ} (hh' : h' < h) (hhH : h < H) :
    ∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
      (if n ∈ Finset.Ico (M - h) (N + 1 - h) then e (φ (n + h)) else 0)
        * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - h') (N + 1 - h')
            then e (φ (n + h')) else 0)
      = ∑ m ∈ Finset.Ico M (N + 1 - (h - h')),
          e (φ (m + (h - h')) - φ m) := by
  classical
  have h1 : ∀ n, (if n ∈ Finset.Ico (M - h) (N + 1 - h)
        then e (φ (n + h)) else 0)
      * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - h') (N + 1 - h')
          then e (φ (n + h')) else 0)
      = if n ∈ Finset.Ico (M - h') (N + 1 - h)
          then e (φ (n + h) - φ (n + h')) else 0 := by
    intro n
    by_cases hn1 : n ∈ Finset.Ico (M - h) (N + 1 - h)
    · by_cases hn2 : n ∈ Finset.Ico (M - h') (N + 1 - h')
      · rw [if_pos hn1, if_pos hn2, e_mul_conj, if_pos ?_]
        rw [Finset.mem_Ico] at hn1 hn2 ⊢
        omega
      · rw [if_pos hn1, if_neg hn2, map_zero, mul_zero, if_neg ?_]
        rw [Finset.mem_Ico] at hn1 hn2 ⊢
        omega
    · rw [if_neg hn1, zero_mul, if_neg ?_]
      rw [Finset.mem_Ico] at hn1 ⊢
      omega
  rw [Finset.sum_congr rfl (fun n _ => h1 n), Finset.sum_ite_mem]
  have h2 : Finset.Ico (M - H + 1) (N + 1)
      ∩ Finset.Ico (M - h') (N + 1 - h)
      = Finset.Ico (M - h') (N + 1 - h) := by
    rw [Finset.inter_eq_right]
    intro n hn
    rw [Finset.mem_Ico] at hn ⊢
    omega
  rw [h2]
  have h4 : ∑ m ∈ Finset.Ico M (N + 1 - (h - h')),
      e (φ (m + (h - h')) - φ m)
      = ∑ n ∈ Finset.Ico (M - h') (N + 1 - h),
          e (φ (n + h) - φ (n + h')) := by
    have h5 : Finset.Ico M (N + 1 - (h - h'))
        = (Finset.Ico (M - h') (N + 1 - h)).map
            (addRightEmbedding h') := by
      rw [Finset.map_add_right_Ico]
      congr 1 <;> omega
    rw [h5, Finset.sum_map]
    refine Finset.sum_congr rfl fun n _ => ?_
    have h6 : (addRightEmbedding h') n = n + h' := rfl
    rw [h6]
    have h7 : n + h' + (h - h') = n + h := by omega
    rw [h7]
  exact h4.symm

/-- The diagonal pair sums are bounded by the window length. -/
private theorem weyl_diag {φ : ℕ → ℝ} {M N H : ℕ} (hHM : H ≤ M) (hMN : M ≤ N)
    {h : ℕ} (hhH : h < H) :
    ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
      (if n ∈ Finset.Ico (M - h) (N + 1 - h) then e (φ (n + h)) else 0)
        * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - h) (N + 1 - h)
            then e (φ (n + h)) else 0)‖
      ≤ ((N + 1 - M : ℕ) : ℝ) := by
  classical
  have h1 : ∀ n, (if n ∈ Finset.Ico (M - h) (N + 1 - h)
        then e (φ (n + h)) else 0)
      * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - h) (N + 1 - h)
          then e (φ (n + h)) else 0)
      = if n ∈ Finset.Ico (M - h) (N + 1 - h) then 1 else 0 := by
    intro n
    by_cases hn : n ∈ Finset.Ico (M - h) (N + 1 - h)
    · rw [if_pos hn, if_pos hn, e_mul_conj, sub_self, e_zero]
    · rw [if_neg hn, zero_mul, if_neg hn]
  rw [Finset.sum_congr rfl (fun n _ => h1 n), Finset.sum_ite_mem,
    Finset.sum_const, nsmul_eq_mul, mul_one]
  have h2 : Finset.Ico (M - H + 1) (N + 1)
      ∩ Finset.Ico (M - h) (N + 1 - h)
      = Finset.Ico (M - h) (N + 1 - h) := by
    rw [Finset.inter_eq_right]
    intro n hn
    rw [Finset.mem_Ico] at hn ⊢
    omega
  rw [h2, Nat.card_Ico, Complex.norm_natCast]
  have h3 : N + 1 - h - (M - h) ≤ N + 1 - M := by omega
  exact_mod_cast h3

/-- Swapping the roles in a conjugate-pair sum preserves the norm. -/
private theorem norm_sum_mul_conj_comm (W : Finset ℕ) (A B : ℕ → ℂ) :
    ‖∑ n ∈ W, A n * (starRingEnd ℂ) (B n)‖
      = ‖∑ n ∈ W, B n * (starRingEnd ℂ) (A n)‖ := by
  have h1 : ∑ n ∈ W, B n * (starRingEnd ℂ) (A n)
      = (starRingEnd ℂ) (∑ n ∈ W, A n * (starRingEnd ℂ) (B n)) := by
    rw [map_sum]
    refine Finset.sum_congr rfl fun n _ => ?_
    rw [map_mul, Complex.conj_conj, mul_comm]
  rw [h1, RCLike.norm_conj]

/-- **Weyl differencing** (the van der Corput A-process): the square of an
exponential sum is controlled by the difference-phase sums. -/
private theorem weyl_differencing {φ : ℕ → ℝ} {M N H : ℕ}
    (hH : 1 ≤ H) (hHM : H ≤ M) (hMN : M ≤ N) :
    ‖∑ n ∈ Finset.Ico M (N + 1), e (φ n)‖ ^ 2
      ≤ (((N + 1 - M : ℕ) : ℝ) + H) / H
        * (((N + 1 - M : ℕ) : ℝ)
          + 2 * ∑ g ∈ Finset.Ico 1 H,
              ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖) := by
  classical
  have hHS : (H : ℂ) * ∑ n ∈ Finset.Ico M (N + 1), e (φ n)
      = ∑ n ∈ Finset.Ico (M - H + 1) (N + 1), ∑ h ∈ Finset.range H,
          (if n ∈ Finset.Ico (M - h) (N + 1 - h)
            then e (φ (n + h)) else 0) := by
    rw [Finset.sum_comm]
    have h1 : ∀ h ∈ Finset.range H,
        ∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
          (if n ∈ Finset.Ico (M - h) (N + 1 - h)
            then e (φ (n + h)) else 0)
        = ∑ n ∈ Finset.Ico M (N + 1), e (φ n) := by
      intro h hh
      rw [Finset.mem_range] at hh
      rw [Finset.sum_ite_mem]
      have h2 : Finset.Ico (M - H + 1) (N + 1)
          ∩ Finset.Ico (M - h) (N + 1 - h)
          = Finset.Ico (M - h) (N + 1 - h) := by
        rw [Finset.inter_eq_right]
        intro n hn
        rw [Finset.mem_Ico] at hn ⊢
        omega
      rw [h2, ← sum_shift φ h (by omega) (by omega)]
    rw [Finset.sum_congr rfl h1, Finset.sum_const, Finset.card_range,
      nsmul_eq_mul]
  -- the pair-sum bound
  have hpairs : ∑ p ∈ (Finset.range H) ×ˢ (Finset.range H),
      ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
        (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
          then e (φ (n + p.1)) else 0)
        * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
            then e (φ (n + p.2)) else 0)‖
      ≤ (H : ℝ) * ((N + 1 - M : ℕ) : ℝ)
        + 2 * ((H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
            ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖) := by
    rw [← Finset.sum_filter_add_sum_filter_not
      ((Finset.range H) ×ˢ (Finset.range H)) (fun p => p.2 < p.1),
      ← Finset.sum_filter_add_sum_filter_not
      (((Finset.range H) ×ˢ (Finset.range H)).filter
        (fun p => ¬ p.2 < p.1)) (fun p => p.1 < p.2)]
    have hS1 : ∑ p ∈ ((Finset.range H) ×ˢ (Finset.range H)).filter
        (fun p => p.2 < p.1),
        ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
          (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
            then e (φ (n + p.1)) else 0)
          * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
              then e (φ (n + p.2)) else 0)‖
        ≤ (H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
            ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖ := by
      have he : ∀ p ∈ ((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => p.2 < p.1),
          ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
            (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
              then e (φ (n + p.1)) else 0)
            * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
                then e (φ (n + p.2)) else 0)‖
          = ‖∑ m ∈ Finset.Ico M (N + 1 - (p.1 - p.2)),
              e (φ (m + (p.1 - p.2)) - φ m)‖ := by
        intro p hp
        rw [Finset.mem_filter, Finset.mem_product, Finset.mem_range,
          Finset.mem_range] at hp
        rw [weyl_offdiag hHM hMN hp.2 hp.1.1]
      rw [Finset.sum_congr rfl he]
      exact sum_offdiag_le
        (fun g => ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖)
        (fun g => norm_nonneg _)
    have hS2 : ∑ p ∈ ((((Finset.range H) ×ˢ (Finset.range H)).filter
        (fun p => ¬ p.2 < p.1)).filter (fun p => p.1 < p.2)),
        ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
          (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
            then e (φ (n + p.1)) else 0)
          * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
              then e (φ (n + p.2)) else 0)‖
        ≤ (H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
            ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖ := by
      have he : ∀ p ∈ (((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => ¬ p.2 < p.1)).filter (fun p => p.1 < p.2),
          ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
            (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
              then e (φ (n + p.1)) else 0)
            * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
                then e (φ (n + p.2)) else 0)‖
          = ‖∑ m ∈ Finset.Ico M (N + 1 - (p.2 - p.1)),
              e (φ (m + (p.2 - p.1)) - φ m)‖ := by
        intro p hp
        rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_product,
          Finset.mem_range, Finset.mem_range] at hp
        rw [norm_sum_mul_conj_comm, weyl_offdiag hHM hMN hp.2 hp.1.1.2]
      rw [Finset.sum_congr rfl he]
      have hswap : ∑ p ∈ (((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => ¬ p.2 < p.1)).filter (fun p => p.1 < p.2),
          ‖∑ m ∈ Finset.Ico M (N + 1 - (p.2 - p.1)),
            e (φ (m + (p.2 - p.1)) - φ m)‖
          = ∑ p ∈ ((Finset.range H) ×ˢ (Finset.range H)).filter
            (fun p => p.2 < p.1),
            ‖∑ m ∈ Finset.Ico M (N + 1 - (p.1 - p.2)),
              e (φ (m + (p.1 - p.2)) - φ m)‖ := by
        refine Finset.sum_nbij' Prod.swap Prod.swap ?_ ?_ ?_ ?_ ?_
        · intro p hp
          rw [Finset.mem_filter, Finset.mem_filter,
            Finset.mem_product] at hp
          rw [Finset.mem_filter, Finset.mem_product]
          exact ⟨⟨hp.1.1.2, hp.1.1.1⟩, hp.2⟩
        · intro p hp
          rw [Finset.mem_filter, Finset.mem_product] at hp
          rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_product]
          simp only [Prod.fst_swap, Prod.snd_swap]
          have h21 := hp.2
          exact ⟨⟨⟨hp.1.2, hp.1.1⟩, by omega⟩, hp.2⟩
        · intro p _
          exact Prod.swap_swap p
        · intro p _
          exact Prod.swap_swap p
        · intro p _
          rfl
      rw [hswap]
      exact sum_offdiag_le
        (fun g => ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖)
        (fun g => norm_nonneg _)
    have hS3 : ∑ p ∈ ((((Finset.range H) ×ˢ (Finset.range H)).filter
        (fun p => ¬ p.2 < p.1)).filter (fun p => ¬ p.1 < p.2)),
        ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
          (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
            then e (φ (n + p.1)) else 0)
          * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
              then e (φ (n + p.2)) else 0)‖
        ≤ (H : ℝ) * ((N + 1 - M : ℕ) : ℝ) := by
      have hbound : ∀ p ∈ (((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => ¬ p.2 < p.1)).filter (fun p => ¬ p.1 < p.2),
          ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
            (if n ∈ Finset.Ico (M - p.1) (N + 1 - p.1)
              then e (φ (n + p.1)) else 0)
            * (starRingEnd ℂ) (if n ∈ Finset.Ico (M - p.2) (N + 1 - p.2)
                then e (φ (n + p.2)) else 0)‖
          ≤ ((N + 1 - M : ℕ) : ℝ) := by
        intro p hp
        rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_product,
          Finset.mem_range, Finset.mem_range] at hp
        have hpe : p.2 = p.1 := by omega
        rw [hpe]
        exact weyl_diag hHM hMN hp.1.1.1
      refine le_trans (Finset.sum_le_sum hbound) ?_
      rw [Finset.sum_const, nsmul_eq_mul]
      refine mul_le_mul_of_nonneg_right ?_ (Nat.cast_nonneg _)
      have hcard : ((((Finset.range H) ×ˢ (Finset.range H)).filter
          (fun p => ¬ p.2 < p.1)).filter (fun p => ¬ p.1 < p.2)).card
          ≤ H := by
        have hinj := Finset.card_le_card_of_injOn Prod.fst
          (s := (((Finset.range H) ×ˢ (Finset.range H)).filter
            (fun p => ¬ p.2 < p.1)).filter (fun p => ¬ p.1 < p.2))
          (t := Finset.range H) ?_ ?_
        · rw [Finset.card_range] at hinj
          exact hinj
        · intro p hp
          rw [Finset.mem_coe, Finset.mem_filter, Finset.mem_filter,
            Finset.mem_product] at hp
          rw [Finset.mem_coe]
          exact hp.1.1.1
        · intro p hp q hq hpq
          rw [Finset.mem_coe, Finset.mem_filter, Finset.mem_filter] at hp hq
          have h5 : p.2 = p.1 := by omega
          have h6 : q.2 = q.1 := by omega
          refine Prod.ext hpq ?_
          rw [h5, h6, hpq]
      exact_mod_cast hcard
    linarith [hS1, hS2, hS3]
  -- assemble
  have hcs := norm_sum_sq_le_card_mul (Finset.Ico (M - H + 1) (N + 1))
    (fun n => ∑ h ∈ Finset.range H,
      (if n ∈ Finset.Ico (M - h) (N + 1 - h) then e (φ (n + h)) else 0))
  have hex := sum_norm_sq_expand (Finset.Ico (M - H + 1) (N + 1))
    (Finset.range H)
    (fun n h => if n ∈ Finset.Ico (M - h) (N + 1 - h)
      then e (φ (n + h)) else 0)
  beta_reduce at hcs hex
  have hWcard : ((Finset.Ico (M - H + 1) (N + 1)).card : ℝ)
      = ((N + 1 - M : ℕ) : ℝ) + (H : ℝ) - 1 := by
    rw [Nat.card_Ico]
    have h8 : N + 1 - (M - H + 1) = (N + 1 - M) + H - 1 := by omega
    rw [h8]
    have h9 : (N + 1 - M) + H - 1 + 1 = (N + 1 - M) + H := by omega
    have h10 := congrArg (fun k : ℕ => (k : ℝ)) h9
    push_cast at h10
    linarith
  have hchain : (H : ℝ) ^ 2
      * ‖∑ n ∈ Finset.Ico M (N + 1), e (φ n)‖ ^ 2
      ≤ (((N + 1 - M : ℕ) : ℝ) + (H : ℝ) - 1)
        * ((H : ℝ) * ((N + 1 - M : ℕ) : ℝ)
          + 2 * ((H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
              ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖)) := by
    have h11 : (H : ℝ) ^ 2
        * ‖∑ n ∈ Finset.Ico M (N + 1), e (φ n)‖ ^ 2
        = ‖(H : ℂ) * ∑ n ∈ Finset.Ico M (N + 1), e (φ n)‖ ^ 2 := by
      rw [norm_mul, Complex.norm_natCast, mul_pow]
    rw [h11, hHS]
    calc ‖∑ n ∈ Finset.Ico (M - H + 1) (N + 1), ∑ h ∈ Finset.range H,
          (if n ∈ Finset.Ico (M - h) (N + 1 - h)
            then e (φ (n + h)) else 0)‖ ^ 2
        ≤ ((Finset.Ico (M - H + 1) (N + 1)).card : ℝ)
          * ∑ n ∈ Finset.Ico (M - H + 1) (N + 1),
            ‖∑ h ∈ Finset.range H,
              (if n ∈ Finset.Ico (M - h) (N + 1 - h)
                then e (φ (n + h)) else 0)‖ ^ 2 := hcs
      _ ≤ ((Finset.Ico (M - H + 1) (N + 1)).card : ℝ)
          * ((H : ℝ) * ((N + 1 - M : ℕ) : ℝ)
            + 2 * ((H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
                ‖∑ m ∈ Finset.Ico M (N + 1 - g),
                  e (φ (m + g) - φ m)‖)) := by
          refine mul_le_mul_of_nonneg_left (le_trans hex hpairs)
            (Nat.cast_nonneg _)
      _ = (((N + 1 - M : ℕ) : ℝ) + (H : ℝ) - 1)
          * ((H : ℝ) * ((N + 1 - M : ℕ) : ℝ)
            + 2 * ((H : ℝ) * ∑ g ∈ Finset.Ico 1 H,
                ‖∑ m ∈ Finset.Ico M (N + 1 - g),
                  e (φ (m + g) - φ m)‖)) := by
          rw [hWcard]
  have hHpos : (0 : ℝ) < (H : ℝ) := by exact_mod_cast hH
  have hFpos : (0 : ℝ) ≤ ∑ g ∈ Finset.Ico 1 H,
      ‖∑ m ∈ Finset.Ico M (N + 1 - g), e (φ (m + g) - φ m)‖ :=
    Finset.sum_nonneg fun g _ => norm_nonneg _
  have hLpos : (0 : ℝ) ≤ ((N + 1 - M : ℕ) : ℝ) := Nat.cast_nonneg _
  rw [div_mul_eq_mul_div, le_div_iff₀ hHpos]
  have hs0 : (0 : ℝ) ≤ ‖∑ n ∈ Finset.Ico M (N + 1), e (φ n)‖ ^ 2 := by
    positivity
  nlinarith [hchain, hHpos, hFpos, hLpos, hs0]


private theorem e_eq_fourierChar (x : ℝ) : e x = (𝐞 x : ℂ) := rfl

/-- Van der Corput differencing, with an adjustable number of shifts. The restriction
`H ≤ a` keeps the shifted windows in the natural numbers; no restriction on the
length of the interval is necessary. -/
theorem norm_sum_Icc_sq_le_correlations {φ : ℕ → ℝ} {a b H : ℕ}
    (hH : 1 ≤ H) (hHa : H ≤ a) (hab : a ≤ b) :
    ‖∑ n ∈ Finset.Icc a b, (𝐞 (φ n) : ℂ)‖ ^ 2 ≤
      (((b + 1 - a : ℕ) : ℝ) + H) / H *
        (((b + 1 - a : ℕ) : ℝ) + 2 * ∑ j ∈ Finset.Ico 1 H,
          ‖∑ n ∈ Finset.Ico a (b + 1 - j), (𝐞 (φ (n + j) - φ n) : ℂ)‖) := by
  simpa only [e_eq_fourierChar, Finset.Ico_add_one_right_eq_Icc] using
    weyl_differencing (φ := φ) hH hHa hab

/-! ## Dyadic selection (proved for Expdb) -/

/-- Dyadic pigeonholing costs only the number of logarithmic ranges. Each range
contains at most its lower endpoint's number of shifts. -/
theorem exists_dyadic_weighted_bound (C : ℕ → ℝ) (hC : ∀ j, 0 ≤ C j)
    {H : ℕ} (hH : 1 ≤ H) :
    ∃ j : ℕ, 1 ≤ j ∧ j ≤ H ∧
      (∑ k ∈ Finset.Ico 1 H, C k) ≤ (Nat.log 2 H + 1 : ℕ) * ((j : ℝ) * C j) := by
  classical
  by_cases hH1 : H = 1
  · refine ⟨1, le_rfl, by omega, ?_⟩
    simp only [hH1, Finset.Ico_self, Finset.sum_empty]
    exact mul_nonneg (Nat.cast_nonneg _) (by simpa using hC 1)
  let J := Finset.Ico 1 H
  obtain ⟨j, hj, hmax⟩ := Finset.exists_max_image J (fun k ↦ (k : ℝ) * C k)
    ⟨1, by simp [J]; omega⟩
  let M : ℝ := (j : ℝ) * C j
  have hM : 0 ≤ M := mul_nonneg (Nat.cast_nonneg _) (hC j)
  have hbin : ∀ k : ℕ, (∑ n ∈ J with Nat.log 2 n = k, C n) ≤ M := by
    intro k
    let B := J.filter (fun n ↦ Nat.log 2 n = k)
    have hB : B ⊆ Finset.Ico (2 ^ k) (2 ^ (k + 1)) := by
      intro n hn
      obtain ⟨hnJ, hnlog⟩ := Finset.mem_filter.1 hn
      have hnpos : n ≠ 0 := by have := (Finset.mem_Ico.1 hnJ).1; omega
      rw [Finset.mem_Ico]
      constructor
      · simpa [hnlog] using Nat.pow_log_le_self 2 hnpos
      · simpa [hnlog] using Nat.lt_pow_succ_log_self (by decide : 1 < 2) n
    have hcard : (B.card : ℝ) ≤ (2 : ℝ) ^ k := by
      have hc := Finset.card_le_card hB
      rw [Nat.card_Ico, pow_succ] at hc
      have hc' : B.card ≤ 2 ^ k := by omega
      exact_mod_cast hc'
    have hpow : (0 : ℝ) < 2 ^ k := by positivity
    have hentry : ∀ n ∈ B, C n ≤ M / (2 : ℝ) ^ k := by
      intro n hn
      have hnJ := (Finset.mem_filter.1 hn).1
      have hlower : (2 : ℝ) ^ k ≤ n := by exact_mod_cast (Finset.mem_Ico.1 (hB hn)).1
      apply (le_div_iff₀ hpow).2
      calc
        C n * 2 ^ k ≤ (n : ℝ) * C n := by nlinarith [hC n]
        _ ≤ M := hmax n hnJ
    calc
      _ ≤ ∑ _n ∈ B, M / (2 : ℝ) ^ k := Finset.sum_le_sum hentry
      _ = (B.card : ℝ) * (M / (2 : ℝ) ^ k) := by simp
      _ ≤ (2 : ℝ) ^ k * (M / (2 : ℝ) ^ k) :=
        mul_le_mul_of_nonneg_right hcard (div_nonneg hM hpow.le)
      _ = M := by field_simp
  have hfiber := Finset.sum_fiberwise_of_maps_to
    (s := J) (t := Finset.range (Nat.log 2 H + 1)) (g := Nat.log 2)
    (fun n hn ↦ Finset.mem_range.2 (by
      exact Nat.lt_succ_of_le (Nat.log_mono_right (Finset.mem_Ico.1 hn).2.le))) C
  refine ⟨j, (Finset.mem_Ico.1 hj).1, (Finset.mem_Ico.1 hj).2.le, ?_⟩
  calc
    _ = ∑ k ∈ Finset.range (Nat.log 2 H + 1), ∑ n ∈ J with Nat.log 2 n = k, C n :=
      hfiber.symm
    _ ≤ ∑ _k ∈ Finset.range (Nat.log 2 H + 1), M :=
      Finset.sum_le_sum fun k _ ↦ hbin k
    _ = _ := by simp [M]

end

end Expdb.VanDerCorput
