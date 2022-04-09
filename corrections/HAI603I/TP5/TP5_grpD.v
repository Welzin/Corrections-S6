Require Import Arith.
Require Export List.
Open Scope list_scope.
Import ListNotations.
Require Import Lia.

(**** Exercice 1 : Spécification ****)

(* Question 1 *)
Inductive is_perm : list nat -> list nat ->Prop:=
 | is_perm_cons : forall (l1 l2 :list nat ) (v:nat ), 
  is_perm l1 l2 -> is_perm (v::l1) (v::l2)
 | is_perm_append : forall (l1 :list nat ) (v:nat ) ,
  is_perm (v::l1 ) (l1 ++ v::nil) 
|is_perm_trans : forall (l1 l2 l3 :list nat ) ,
  is_perm l1 l2 -> is_perm l2 l3 -> is_perm l1 l3
|is_perm_refl : forall  (l1 : list nat ), is_perm l1 l1
|is_perm_sym : forall  (l1 l2 : list nat ), is_perm l1 l2 ->
is_perm l2 l1 .

(* Question  2 *)
Lemma is_perm_ex1 : is_perm [1;2;3] [3;2;1].
Proof.
apply (is_perm_trans [1; 2; 3] ([2; 3] ++ [1]) [3; 2; 1]).
 (* is_perm [1; 2; 3] ([2; 3] ++ [1]) *)
  -simpl.
  apply is_perm_append .
(* is_perm ([2; 3] ++ [1]) [3; 2; 1] *)
  -apply (is_perm_trans ([2;3]++[1]) ( [3]++[1]++[2]) [3; 2; 1]).
   (*is_perm ([2; 3] ++ [1]) ([3] ++ [1] ++ [2]) *)
    *simpl.
    apply is_perm_append.
   (*is_perm ([3] ++ [1] ++ [2]) [3; 2; 1]*)
    *simpl.
    apply is_perm_cons.
    apply is_perm_append.
Qed.

(* Question  3 *)
Inductive is_sorted : list nat ->Prop :=
|is_sorted_O : is_sorted nil 
|is_sorted_sing :forall (v1:nat),  is_sorted [v1]
|is_sorted_cons : forall (l1 : list nat ) (v1 v2 : nat ) ,
 v1 <= v2 -> is_sorted (v2::l1)->is_sorted(v1::v2::l1) .

(* Question 4 *)
Lemma is_sorted_ex1 : is_sorted [1; 2; 3].
Proof.
apply is_sorted_cons.
(* 1 < = 2 *) 
apply Nat.le_succ_diag_r.
(* is_sorted [2; 3] *) 
apply is_sorted_cons.
(* 2 < = 3 *)
apply Nat.le_succ_diag_r.
(* is_sorted [3] *) 
apply is_sorted_sing.
Qed.


Ltac is_sorted_tac :=
(* Cas de base *)
try apply is_sorted_O; 
(* Cas reccurent *)
repeat ( apply is_sorted_cons || apply Nat.le_succ_diag_r ); (*decomposition liste avec :  _ <= _ *)
try apply is_sorted_sing. (* ne reste plus que le singleton *) 


Lemma is_sorted_testTact_0 : is_sorted [].
Proof.
is_sorted_tac.
Qed.

Lemma is_sorted_testTact_1 : is_sorted [1].
Proof.
is_sorted_tac.
Qed.

Lemma is_sorted_testTact_n : is_sorted [1; 2; 3; 4; 5].
Proof.
is_sorted_tac.
Qed.




(**** Exercice 2 : Fonction de tri par insertion ****)
(* Question 1 *)

Check le_dec.
Print sumbool.

(* ce lemme peut être vu comme une fonction qui rend si n ≤ m left sinon right  *)

Definition le_10 ( n : nat ) : bool :=
match ( le_dec n 10 ) with
| left _ => true
| right _ => false
end .

Eval compute in (le_10 5). (**true**)
Eval compute in (le_10 15). (** false **)

(* Question 2 *)
Fixpoint insert (x : nat) (l : list nat) {struct l} : list nat :=
match l with
|nil => [x]
|h::t => 
  match (le_dec x h) with 
    |left _ => x::h::t 
    | right _ => h::(insert x t )
  end
end.

(* Question 3 *)
Eval compute in (insert 3 [1;2;3;5]). (** [1;2;3;3;5] **)

Lemma insert_ex1 : insert 3 [1;2;3;5] = [1;2;3;3;5].
Proof.
simpl.
reflexivity.
Qed.


(* Question 4 *)
Fixpoint isort (l : list nat) : list nat :=
match l with 
|nil => nil
| v::l => insert v (isort l)
end.

(* Question 5 *)
Eval compute in (isort [5; 4; 3; 2; 1]).



(**** Exercice 3 : Preuve de correction ****)

Lemma isort_ex1 : isort [5; 4; 3; 2; 1] = [1; 2; 3; 4; 5].
Proof.
simpl.
reflexivity.
Qed.

(* Question 1 *)
Lemma head_is_perm : forall (x1 x2 : nat) (l : list nat),
is_perm (x1 :: x2 :: l) (x2 :: x1 :: l).
Proof.
intros x1 x2 l.
apply (is_perm_trans (x1 :: x2 :: l)  ((x2 :: l )++[x1]) (x2 :: x1 :: l)).
 (* Cas is_perm (x1 :: x2 :: l) ((x2 :: l) ++ [x1]) *)
  -apply is_perm_append.
 (* Cas is_perm ((x2 :: l) ++ [x1]) (x2 :: x1 :: l) *)
  -simpl.
  apply is_perm_cons.
  apply is_perm_sym.
  apply is_perm_append.
Qed.

(* Question 2 *)
Lemma insert_is_perm : forall (x : nat) (l : list nat),
is_perm (x::l) (insert x l).
Proof.
intros x l .
induction l as [| a l].
 (* Cas de base *)
  - simpl.
  apply is_perm_refl.
 (* Cas inductif *)
  -simpl.
  destruct (le_dec x a)  as [ x_inf_a | x_ninf_a ].
   (* Cas x <= a *)
    *apply is_perm_refl.
   (* Cas ~ x <= a *)
    *apply ( is_perm_trans (x :: a :: l) (a::x::l) (a :: insert x l)).
      (* Cas is_perm (x :: a :: l) (a :: x :: l) *)
        ** apply head_is_perm.
      (* Cas is_perm (a :: x :: l) (a :: insert x l) *)
        **apply is_perm_cons;assumption.
Qed.

(* Question 3 *)
Lemma insert_is_sorted : forall (x:nat ) (l:list nat ) ,
is_sorted l ->is_sorted (insert x l).
Proof.
intros x l H.
induction H as [ |  v1 | l0 v1 v2  v1_inf_v2 isS_v2l0 IHisS ]. (* on fait une induction sur une relation donc on doit forunir les noms des variables pour chacun des constructeur de la relation H *)
 (* Cas is_sorted_nil *)
  -simpl.
  apply is_sorted_sing.
 (* Cas is_sorted_sing *)
  -simpl.
  destruct (le_dec x v1) as [ x_inf_v1 | x_ninf_v1 ].
   (* Cas  x <= v1 *)
    *apply is_sorted_cons.
    apply x_inf_v1.
    apply is_sorted_sing.
  (* Cas ~ x <= v1*)
    *apply is_sorted_cons.
    generalize x_ninf_v1.
    rewrite -> Nat.nle_gt. 
    apply Nat.lt_le_incl.
    apply is_sorted_sing.
 (* Cas is_sorted_cons*)
  -simpl.
  destruct (le_dec x v1) as [ x_inf_v1 | x_ninf_v1 ].
   (* Cas  x <= v1 *)
    *apply is_sorted_cons.
    assumption.
    apply is_sorted_cons.
    assumption.
    apply isS_v2l0 .  (* <=> assumption *)
   (* Cas  ~ x <= v1 *)
    *rewrite -> Nat.nle_gt in x_ninf_v1 ;apply Nat.lt_le_incl in x_ninf_v1.
    destruct (le_dec x v2) as [ x_inf_v2 | x_ninf_v2 ] eqn:lol1 .
     (* Cas   v1 <= x et x <= v2 *)
      **apply is_sorted_cons.
      assumption.
      apply is_sorted_cons.
      assumption.
      apply isS_v2l0 .
     (* Cas   v1 <= x et ~ x <= v2 *)
      **apply is_sorted_cons.
      assumption.
      simpl in IHisS. 
      destruct (le_dec x v2) as [ x_inf_v2' | x_ninf_v2' ] eqn:lol2.
       (* Cas v1 <= x et ~ x <= v2 et  x <= v2 *)
        ***discriminate. (* lia *) (*voir note *)
       (* Cas v1 <= x et ~ x <= v2 et  ~x <= v2 *)
        ***assumption.
Qed.


(* Question 4 *)
Lemma isort_correct : forall ( l l' : list nat ) ,
l' = isort l -> is_perm l l' /\  is_sorted l' .
Proof.
induction l as [| a l1 IHl].
(*Cas de base *)
intros l' H0;subst;simpl.
split.
(* Cas is_perm *)
apply is_perm_refl.
(* Cas is_sorted *)
apply is_sorted_O.
(*Cas inductif *)
intros l' H1;subst;simpl.
split.
(* Cas is_perm  *) 
apply (is_perm_trans (a :: l1) (a ::isort l1)  (insert a (isort l1))).
apply is_perm_cons.
apply IHl.
trivial.
apply insert_is_perm.
(* Cas is_sorted  *) 
apply insert_is_sorted.
apply IHl.
trivial.
Qed.

(**** Exercice 4 : Preuve de complétude ****)

(* Question 1 *)

Lemma is_sorted_cons_inv : forall (l : list nat) (a : nat), is_sorted (a :: l) -> is_sorted l.
Proof.
intros l  a .
induction l as [ | a0 l  IHl].
intro ;apply is_sorted_O.
intro H.
inversion H as [ |  |lc v1 v2  v1_inf_v2 isS_v2lc  ]. (* voir constructeur v1_inf_v2->isS_v2l-> is_sorted (v1 :: v2 :: lc) *)
apply isS_v2lc. (*  H0 : v1 = a ,H1 : v2 = a0  + subst donc on peut utiliser isS_v2l à la place de isS_a0l *) 
Qed.




(* Question 2 *)

Lemma is_sorted_app :  forall (l : list nat) (a : nat), is_sorted (l ++ [a]) -> is_sorted l.
Proof.
induction l as [ | a0 l1  IHl1].
(*Cas de base *)
intros.
apply is_sorted_O.
(*Cas Inductif *)
intros a H.
induction l1 as [ | a1 l2  IHl2].
(*Cas de base *)
apply is_sorted_sing.
(*Cas Inductif *)
apply is_sorted_cons.
inversion H as [ |  |lc v1 v2  v1_inf_v2 isS_v2lc  ].
assumption.
apply (IHl1 a) .
(*inversion H.
assumption.*)
rewrite <-app_comm_cons in H.
apply (is_sorted_cons_inv  ((a1 :: l2) ++ [a]) a0 )in H.
assumption.
Qed.

(* Correction groupe D *)

Lemma is_sorted_app2 :
forall l a, is_sorted (l ++ [a]) -> is_sorted l.
Proof.
induction l.
(* Cas de base *)
- intros; apply is_sorted_O.
(* Cas inductif *)
- intros x iss_app.
(* On transforme is_sorted ((a :: l) ++ [x]) en
is_sorted (a :: l ++ [x]) avant de raisonner sur la definition
de is_sorted en utilisant la tactique inversion. *)
Check app_comm_cons.
rewrite <- app_comm_cons in iss_app.
inversion iss_app.
(* Cas is_sorted_sing *)
(* Contradiction par [] = l ++ [x] *)
+ elimtype False.
generalize H1.
apply app_cons_not_nil.
(* Cas is_sorted_cons. *)
+ (* Il faut raisonner sur la forme de l. *)
destruct l.
(* Cas nil *)
-- apply is_sorted_sing.
(* Cas cons *)
-- apply is_sorted_cons.
cut (n = v2).
intros eq_nm; rewrite eq_nm; assumption.
rewrite <- app_comm_cons in H0.
injection H0; intros _ eq_mn; rewrite eq_mn; reflexivity.
apply IHl with (a := x).
rewrite <- H0; assumption.
Qed.
(* Fin correction groupe D *)


(* Question 3 *)

Lemma is_sorted_id : forall l, is_sorted l -> isort l = l.
Proof.
induction l as [ | a0 l1  IHl1];simpl.
(*Cas de base *)
intro H0.
trivial.
(*Cas Inductif *)
intro H.
rewrite IHl1.

- induction l1 as [ | a1 l2  IHl2];simpl.
(*Cas de base *)
trivial.
(*Cas Inductif *)
destruct (le_dec a0 a1) as [a0_inf_a1 | a0_ninf_a1] .
(*Cas a0 <= a1 *)
trivial.
(*Cas ~ a0 <= a1 *)
elimtype False.
apply a0_ninf_a1.
(* lia *)
inversion H as [ |  |lc v1 v2  v1_inf_v2 isS_v2lc  ].
apply v1_inf_v2. (* assumption  *)

- apply is_sorted_cons_inv in H.
assumption.
(*
inversion H . (* si is_sorted l1 a été déduit de is_sorted (a0 :: l1) alors forcement par singleton donc l1 =[]   *)
apply is_sorted_O. (* is sorted l1 vue que l1 = [] verifier par _0 *)
apply is_sorted_cons_inv in H.
assumption.*)
Qed.


(* Question 4 *)
Lemma is_sorted_insert_id : forall l a,is_sorted (a::l) -> insert a l = a::l .
Proof.
induction l as [ | a0 l1  IHl1];simpl.
(*Cas de base *)
trivial.
(*Cas Inductif *)
intros a1 H.
destruct (le_dec a1 a0) as [a1_inf_a0 | a1_ninf_a0] .
(*Cas  a1 <= a0 *)
trivial.
(*Cas  ~ a1 <= a0 *)
elimtype False.
apply a1_ninf_a0.
inversion H as [ |  |lc v1 v2  v1_inf_v2 isS_v2lc  ].
apply v1_inf_v2 .
Qed.


(* Question 5 *)

Lemma is_sorted_not_lt_tl_hd : forall l a x , is_sorted (a::l ++[x]) -> x < a ->False.
Proof.
induction l as [ | a0 l1  IHl1].
(* Cas de base *)
intros a x H.
simpl in H.
inversion H as [ |  |lc v1 v2  v1_inf_v2 is_Sv2lc  ]. 
Search (  not (_ > _) ).
apply (le_not_gt a x) in v1_inf_v2.
assumption.
(* Cas Inductif *)
intros.
inversion H as [ |  |lc v1 v2  v1_inf_v2 is_Sv2lc  ]. 
apply (IHl1  a0 x) .
assumption.
Search( _ < _ -> _ <= _ -> _ < _ ).
apply (Nat.lt_le_trans x a a0) .
assumption.
assumption.
Qed.

(* Question 6 *)



Lemma is_sorted_app_insert : forall l a , is_sorted(l++[a]) -> insert a l = l++[a].
Proof.
induction l as [ | a0 l1  IHl1].
(* Cas de base *)
intros.
simpl.
trivial.
(* Cas Inductif *)
intros a H.
simpl.
destruct (le_dec a a0 ) as [ a_inf_a0 | a_ninf_a0].
(* Cas a <= a0 *)
-Search( _ <= _ -> _ < _ \/ _ = _).
apply le_lt_or_eq in a_inf_a0.
inversion a_inf_a0 as [ a_sinf_a0|a_eq_a0].
(* Cas a < a0 *)
+elimtype False.
apply (is_sorted_not_lt_tl_hd l1 a0 a).
rewrite <- app_comm_cons in H.
apply H.
apply a_sinf_a0.
(* Cas a = a0 *)
+rewrite a_eq_a0 .
rewrite <- IHl1.
rewrite is_sorted_insert_id.
trivial.
apply is_sorted_app  in H.
assumption.
rewrite <- app_comm_cons in H.
apply is_sorted_cons_inv  in H.
rewrite <- a_eq_a0 .
assumption.
(* Cas ~ a <= a0 *)
-rewrite IHl1.
trivial.
rewrite <- app_comm_cons in H.
apply is_sorted_cons_inv in H.
assumption.
Qed.

(* Question 7 : *)

Hypothesis insert_comm_assoc : forall l a a0, insert a (insert a0 l) = insert a0 (insert a l) .

Lemma isort_eq_cons_app : forall l a , isort(a::l) =isort(l++[a]).

Proof.
induction l as [ | a1 l1  IHl1].
(* Cas de base *)
intro a0;simpl;trivial.
(* Cas inductif *)
intro a0 ;simpl.
rewrite <- IHl1.
simpl.
apply insert_comm_assoc.
Qed.

(* BONUS : *)
 
Lemma insert_comm_assoc_bonus : forall l a a0, insert a (insert a0 l) = insert a0 (insert a l) .
Proof.
induction l;intros a0 a1;simpl;repeat (destruct le_dec;simpl).
all: try lia. (* cas résolution 1 *)
all: try (cut (a0 = a1); [intros e; subst; reflexivity | lia ]). (* cas résolution 2 *)
all: try reflexivity. (* cas résolution 3 *)
(* cas relation rec  *)
rewrite IHl.
trivial.
Qed.

(* Correction Groupe D : *)

Lemma insert_comm_assoc_bonus2 :
forall l x y, insert x (insert y l) = insert y (insert x l).
Proof.
induction l.
(* Cas de base *)
intros x y;simpl;  destruct le_dec as [ x_inf_y | x_ninf_y ] ; destruct le_dec  as [ y_inf_x | y_ninf_x ];
(lia || reflexivity || apply (Nat.le_antisymm x y) in x_inf_y;subst;trivial).
(* Cas inductif *)
intros x y; simpl.
repeat (simpl;destruct le_dec);(* voir fin du TP pour comprendre nb buts généré *)
( lia || reflexivity || 
(cut (x = y);intros;subst;trivial;apply Nat.le_antisymm ;assumption) (* voir fin TP pour erreur interessante *) 
||(rewrite IHl; reflexivity )).
(* TODO : pourquoi paranthéser  les exp est nécessaire ?  *)
Qed.


(* Question 8 : *)

Lemma is_perm_eq_isort : forall l m , is_perm l m -> isort l = isort m.
Proof.
intros l m H.
induction H.
(* Cas cons *)
simpl.
rewrite IHis_perm.
trivial.
(* Cas app *)
apply isort_eq_cons_app.
(* Cas  trans *)
rewrite <- IHis_perm2.
assumption.
(* Cas  refl *)
trivial.
(* Cas sym *)
rewrite IHis_perm.
trivial.
Qed.


(* Question 9 : *)

Lemma is_sort_compl : forall (l l' : list nat),
    is_perm l l' -> is_sorted l' -> isort l = l'.
Proof.
intros l1 l2 H.
induction H. 

    (* Cas cons *)
    intro .
    simpl.
    rewrite IHis_perm.
    apply is_sorted_insert_id in H0.
    assumption.
    apply  is_sorted_cons_inv  in H0.
    assumption.
  
    (*Cas app *)
    intro.
    simpl.
    assert( H' := H).
    apply is_sorted_app_insert in H'.
    apply is_sorted_app in H .
    apply is_sorted_id in H.
    rewrite H.
    rewrite H'.
    trivial.

    (* Cas trans *)
    intro.
    rewrite <- is_perm_eq_isort with (l := l3).
    apply is_sorted_id.
    assumption.
    apply (is_perm_trans l3 l2 l1)  .
    apply is_perm_sym;assumption.
    apply is_perm_sym;assumption.

    
    (* Cas Refl *)
    apply is_sorted_id. (* intro ;apply is_sorted_id ;assumption .*)
    
    (*Cas sym *)
    intro .
    rewrite <- is_sorted_id with( l := l1).
    apply is_perm_eq_isort.
    apply is_perm_sym;assumption.
    assumption. 
Qed.

(* Correction du Groupe D pimpé (voir le moodle pour la preuve  vanilla)   : *)

Lemma is_sort_compl2 : forall (l l' : list nat),
is_perm l l' -> is_sorted l' -> isort l = l'.
Proof.

induction 1 as [  l0 l1 a  |   l0 a  |   l0 l1  l2  is_p_l0l1  |  l0  |   l0 l1 is_p_l0l1 isS   ].
(* Cas cons *)
- simpl.
intros is_s_al1.
inversion is_s_al1 as [ |  |lc v1 v2  v1_inf_v2 is_Sv2lc  ] ; subst.
(* Cas is_sorted_nil, resolu automatiquement par inversion. *)
(* Cas is_sorted_sing *)
rewrite (IHis_perm is_sorted_O).
simpl; reflexivity.
(* Cas is_sorted_cons *)
rewrite IHis_perm.
simpl. elim (le_dec a v2).
reflexivity.
intros nle_am; elimtype False; apply nle_am;apply v1_inf_v2.
apply is_Sv2lc.
(* Cas app *)
- simpl.
intros is_s_l0a.
rewrite is_sorted_id.
rewrite is_sorted_app_insert.
reflexivity.
apply is_s_l0a.
apply ( is_sorted_app l0 a  ). (* apply is_sorted_app with (a := a). *)
apply is_s_l0a.

(* Cas trans *)
- intros isSl2 ;transitivity (isort l1).
apply is_perm_eq_isort.
apply is_p_l0l1.
apply IHis_perm1;apply isSl2.

(* Cas refl *)
- apply is_sorted_id.

(* Cas sym *)
- intros is_s_l0;rewrite <- (is_perm_eq_isort l0 l1) .
apply is_sorted_id;apply is_s_l0 .
apply is_p_l0l1.
Qed.

(* fin correction Groupe D *)

(**** Annexe : ****)

(* Resultat de : "Show Proof." Question 7  Bonus  

   (fun (a : nat) (l0 : list nat)
      (IHl : forall x y : nat, insert x (insert y l0) = insert y (insert x l0))
      (x y : nat) =>
    let s := le_dec y a in
    | left l1 =>
        let s0 := le_dec x y in
        | left l2 =>
            let s1 := le_dec x a in
            | left l3 =>
                let s2 := le_dec y x in
                | left l4 => ?Goal@{l:=l0; l0:=l1; l1:=l2; l2:=l3; l3:=l4} // profondeur : 4 le_dec 
                | right n =>
                    let s3 := le_dec y a in
                    | left l4 => ?Goal0@{l:=l0; l0:=l1; l1:=l2; l2:=l3; l3:=l4} // 5 le_dec
                    | right n0 => ?Goal1@{l:=l0; l0:=l1; l1:=l2; l2:=l3}// 5 le_dec
                    end
                end
            | right n =>
                let s2 := le_dec y a in
                | left l3 => ?Goal2@{l:=l0; l0:=l1; l1:=l2; l2:=l3}// 4 le_dec 
                | right n0 => ?Goal3@{l:=l0; l0:=l1; l1:=l2}// 4 le_dec 
                end
            end
        | right n =>
            let s1 := le_dec x a in
            | left l2 =>
                let s2 := le_dec y x in
                | left l3 => ?Goal4@{l:=l0; l0:=l1; l1:=l2; l2:=l3}// 4 le_dec 
                | right n0 =>
                    let s3 := le_dec y a in
                    | left l3 => ?Goal5@{l:=l0; l0:=l1; l1:=l2; l2:=l3}// 5 le_dec 
                    | right n1 => ?Goal6@{l:=l0; l0:=l1; l1:=l2}// 5 le_dec 
                    end
                end
            | right n0 =>
                let s2 := le_dec y a in
                | left l2 => ?Goal7@{l:=l0; l0:=l1; l1:=l2}// 4 le_dec 
                | right n1 => ?Goal8@{l:=l0; l0:=l1}// 4 le_dec 
                end
            end
        end
    | right n =>
        let s0 := le_dec x a in
        | left l1 =>
            let s1 := le_dec y x in
            | left l2 => ?Goal9@{l:=l0; l0:=l1; l1:=l2}// 3 le_dec  
            | right n0 =>
                let s2 := le_dec y a in
                | left l2 => ?Goal10@{l:=l0; l0:=l1; l1:=l2}// 4 le_dec 
                | right n1 => ?Goal11@{l:=l0; l0:=l1}// 4 le_dec 
                end
            end
        | right n0 =>
            let s1 := le_dec y a in
            | left l1 => ?Goal12@{l:=l0; l0:=l1}// 3 le_dec  
            | right n1 => ?Goal13@{l:=l0}// 3 le_dec  
            end
        end
    end) l)

*)

(* Si on veut essayer avec :
repeat (simpl;destruct le_dec);try( lia || reflexivity || (rewrite IHl; reflexivity ) ||
(cut (x = y);intros;subst;trivial;apply Nat.le_antisymm ;auto)
Alors on doit respecter la precedence décrite ci dessous pour évtier les conflicts :
(*
(* Cas Inductif : Preuve Manuelle *)
cut (x = y);intros  ;subst ;trivial;apply Nat.le_antisymm ;assumption.
reflexivity. (*et cut/auto fonctionne donc precedence (refl ) > precedence (cut (x = y);[...];auto ) *)
lia.
lia.
lia. (* et rw IHl fonctionne donc precedence (lia ) > precedence (rw IHl) *)
reflexivity. (* et cut/auto  *)
lia. (* et fonctionne donc precedence (lia ) > precedence (cut (x = y);[...];auto )  *)
lia.
reflexivity. (* et cut/auto  *)
lia. (* et rw IHl *)
lia.
lia.
reflexivity. 
lia. (* et rw IHl *)
rewrite IHl;reflexivity. (* et cut/auto  fonctionne donc precedence (rw IHl ) > precedence (cut (x = y);[...];auto )  *)
*)

Qed.*)

