use std::marker::PhantomData;

use rand::Rng;

use elliptic_curve_algorithms::field::*;

#[macro_use] use elliptic_curve_algorithms;
use elliptic_curve_algorithms::finite_fields::{Fp, IntegerAsType};
use elliptic_curve_algorithms::finite_fields::integer_mpz::*;

use elliptic_curve_algorithms::elliptic_curves::*;
use elliptic_curve_algorithms::elliptic_curves::fp_elliptic_curves::*;

use elliptic_curve_algorithms::finite_fields::FiniteField;

pub type Integer = gmp::mpz::Mpz;

type PrimeList = Vec<Integer>;

pub struct CSIDHInstance<K>{
    data : PhantomData<K>,
    pub n_primes : usize,
    pub l: PrimeList,
    pub p: Integer,
}


impl<N> CSIDHInstance<Fp<N, Integer>>
where N : IntegerAsType<Integer>{

    fn compute_rhs(x : &Fp<N, Integer>, a : &Fp<N, Integer>) -> Fp<N, Integer>{
        let x_sq = x*x;
        let y = &x_sq*x  + a*&x_sq;
        &y + x
    }

    pub fn new(n_primes : usize, p: Integer, l: PrimeList) -> CSIDHInstance<Fp<N, Integer>>{
        CSIDHInstance{
            data: PhantomData,
            n_primes,
            l, p
        }
    }

    pub fn check_well_defined(self : &CSIDHInstance<Fp<N, Integer>>){
        let L = &self.l;
        let P = &self.p;
        let N_PRIMES = self.n_primes;

        assert_eq!(N_PRIMES, self.l.len());

        let mut prod = Integer::from(1);
        for i in 0..N_PRIMES{
            prod *= L[i].clone();
        }
        assert_eq!(*P, Fp::<N, Integer>::cardinal());
        assert_eq!(*P, Integer::from(4)*prod-Integer::from(1));
    }

    fn is_supersingular(self : &CSIDHInstance<Fp<N, Integer>>, ell : &EllipticCurve<Fp<N, Integer>>) -> bool{
        let L = &self.l;
        let P = &self.p;
        let N_PRIMES = self.n_primes;

        loop{
            let p = EllipticCurve::unsigne_point(ell.sample_point()); // We test multiple point if necessary

            let mut d : Integer = Integer::from(1);

            for i in 0..N_PRIMES{
                let li = L[i].clone();
                let qi = ell.scalar_mult_unsigned((P+1)/li.clone(), p.clone());

                if ell.scalar_mult_unsigned(li.clone(), qi.clone()) != UnsignedProjPoint::infinite_point() {
                    return false;
                }
                if qi != UnsignedProjPoint::infinite_point(){
                    d *= li;
                }
                if &d*&d > &Integer::from(16)*P{
                    return true;
                }
            }
        }
    }

    fn order_naive(ell : &EllipticCurve<Fp<N, Integer>>, p : &UnsignedProjPoint<Fp<N, Integer>>) -> Integer{
        assert!(ell.is_montgomery());

        if p == &UnsignedProjPoint::infinite_point(){
            return Integer::from(1);
        }

        let p_normalized = p.clone().normalize();

        let mut order = Integer::from(2);

        // we start at i = 2 because of special doubling case
        let mut t = ell.x_dbl(p_normalized.clone()); // t = [i]p will iterate over elements of <p>
        let mut t_minus_1 = p_normalized;

        while t != UnsignedProjPoint::infinite_point(){
            order += Integer::from(1);

            let _temp = t.clone();
            t = ell.x_add(t, p.clone(), t_minus_1);
            t_minus_1 = _temp;
        }

        order
    }

    fn isogeny_velu(ell : &EllipticCurve<Fp<N, Integer>>, point : &UnsignedProjPoint<Fp<N, Integer>>, mut q : UnsignedProjPoint<Fp<N, Integer>>, k : Integer) -> Result<(EllipticCurve<Fp<N, Integer>>, UnsignedProjPoint<Fp<N, Integer>>), ()>{
        assert!(ell.is_montgomery());
        assert!(k>=Integer::from(3));
        assert!(&k%2 != Integer::from(0));

        let p = point;

        let mut i = Integer::from(1);
        
        let mut t = p.clone();  // t = [i]p will iterate over elements of <p>
        let mut t_minus_1 = UnsignedProjPoint::infinite_point();

        let mut pi = UnsignedProjPoint::finite_point(Fp::from_int(1));
        let mut sigma = Fp::from_int(0);
        
        let mut projection_x = Fp::from_int(1);
        let mut projection_z = Fp::from_int(1);

        while &Integer::from(2)*&i < k{
            if t.x == Fp::from_int(0){ // point of order 2, we abort
                return Err(());
            }

            pi.x *= t.x.clone();
            pi.z *= t.z.clone();

            sigma += (&t.x*&t.x - &t.z*&t.z)/(&t.z*&t.x);

            projection_x *= &t.x*&q.x - &t.z*&q.z;
            projection_z *= &q.x*&t.z - &t.x*&q.z;

            if i == Integer::from(1){
                let _temp = t.clone();
                t = ell.x_dbl(p.clone());
                t_minus_1 = _temp;
            }else{
                let _temp = t.clone();
                t = ell.x_add(t, p.clone(), t_minus_1);
                t_minus_1 = _temp;
            }
            i += Integer::from(1);
        }

        pi.x *= pi.x.clone();
        pi.z *= pi.z.clone();

        projection_x *= projection_x.clone();
        projection_z *= projection_z.clone();

        sigma *= Fp::from_int(2);

        pi = pi.normalize();
        Ok((
            EllipticCurve::new_montgomery(pi.x*(&ell.a_2 - &(Fp::from_int(3)*sigma))),
            UnsignedProjPoint{
                    x: q.x*projection_x,
                    z: q.z*projection_z
                }
        ))

    }

    pub fn verify_public_key(self : &CSIDHInstance<Fp<N, Integer>>, pk : Fp<N, Integer>) -> bool{
        self.is_supersingular(&EllipticCurve::new_montgomery(pk))
    }

    pub fn naive_class_group_action(self : &CSIDHInstance<Fp<N, Integer>>, pk : Fp<N, Integer>, sk : Vec<i32>) -> Fp<N, Integer>{
        let L = &self.l;
        let P = &self.p;
        let N_PRIMES = &self.n_primes;

        let mut ell = EllipticCurve::new_montgomery(pk);

        let p_plus_1 = P+Integer::from(1);

        for i in 0..sk.len(){
            if sk[i] == 0{
                continue;
            }
            let s = if sk[i]>0{ 1 } else { -1 };

            for _j in 0..(s*sk[i]){
                let mut x = Fp::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));
                let mut p_point = UnsignedProjPoint::finite_point(x.clone());

                let mut q_point = ell.scalar_mult_unsigned(p_plus_1.clone()/L[i].clone(), p_point.clone());

                while CSIDHInstance::compute_rhs(&x, &ell.a_2).legendre_symbol() != s as i8 || q_point == UnsignedProjPoint::infinite_point() {
                    x = Fp::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));

                    p_point = UnsignedProjPoint::finite_point(x.clone());
                    q_point = ell.scalar_mult_unsigned(p_plus_1.clone()/L[i].clone(), p_point);
                }

                let (ell_, _) = Self::isogeny_velu(&ell, &q_point, q_point.clone(), L[i].clone()).unwrap();
                ell = ell_;
            }

        }
        ell.a_2
    }

    pub fn class_group_action(self : &CSIDHInstance<Fp<N, Integer>>, pk : Fp<N, Integer>, mut sk : Vec<i32>) -> Fp<N, Integer>{
        let L = &self.l;
        let P = &self.p;
        let N_PRIMES = &self.n_primes;
        let mut finished_total : [bool; 2] = [false, false];

        let mut k_sign : [Integer; 2]= [Integer::from(1), Integer::from(1)]; // 1 -> >= 0;  0-> <= 0
        let mut s_sign : [Vec<Integer>; 2] = [vec!(), vec!()];
        let mut e_sign : [Vec<i32>; 2] = [vec!(), vec!()];
        let mut finished_sign : [Vec<bool>; 2] = [vec!(), vec!()];

        for i in 0..sk.len(){
            if sk[i] == 0{
            }else if sk[i] > 0{
                e_sign[1].push(sk[i]);
                s_sign[1].push(L[i].clone());
                finished_sign[1].push(false);

                k_sign[1] *= L[i].clone();
            }else{
                e_sign[0].push(sk[i]);
                s_sign[0].push(L[i].clone());
                finished_sign[0].push(false);

                k_sign[0] *= L[i].clone();
            }
        }

        let mut ell = EllipticCurve::new_montgomery(pk);

        while !finished_total[0] || !finished_total[1] {
            // Sample the point
            let x = Fp::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));

            let s = CSIDHInstance::compute_rhs(&x, &ell.a_2).legendre_symbol() as i32;
            if s == 0{
                continue;
            }
            let sign_index = ((s + 1)/2) as usize;
            
            if finished_total[sign_index]{
                continue;
            }
            
            finished_total[sign_index] = true;

            let uns_p = UnsignedProjPoint::finite_point(x);

            let mut k = k_sign[sign_index].clone();

            let mut q_point = ell.scalar_mult_unsigned((P.clone()+Integer::from(1))/k.clone(), uns_p);

            for j in 0..s_sign[sign_index].len(){
                
                // Trick from the original implementation, they started by the big primes
                let i = s_sign[sign_index].len()-1-j;

                if finished_sign[sign_index][i]{
                    continue;
                }
                finished_total[sign_index] = false;

                let li : Integer = s_sign[sign_index][i].clone();

                let r_point = ell.scalar_mult_unsigned(&k/(&li), q_point.clone());

                if &r_point == &UnsignedProjPoint::infinite_point(){
                    // The point sampled had too low l-index
                    continue;
                }

                let (ell_, q_point_) = match Self::isogeny_velu(&ell, &r_point, q_point.clone(), li.clone()){
                    Err(()) => {
                        // This should never happend
                        println!("Erreur");
                        continue;
                        },
                    Ok(pair) => pair
                };

                ell = ell_;
                q_point = q_point_;

                e_sign[sign_index][i] -= s;

                if e_sign[sign_index][i] == 0{
                    finished_sign[sign_index][i] = true;
                }

                k /= li;
            }
        }
        ell.a_2
    }

    pub fn sample_keys(self : &CSIDHInstance<Fp<N, Integer>>, m : i32) -> (Fp<N, Integer>, Vec<i32>){
        let mut rng = rand::thread_rng();
        let N_PRIMES = self.n_primes;
        let mut sk : Vec<i32> = vec!();
        for _i in 0..N_PRIMES{
            sk.push(rng.gen_range(-m, m) as i32);
        }

        let pk = self.class_group_action(Fp::from_int(0), sk.clone());
        (pk, sk)
    }

    pub fn sample_keys_modif(self : &CSIDHInstance<Fp<N, Integer>>, bound : i32) -> (Fp<N, Integer>, Vec<i32>){
        let mut rng = rand::thread_rng();
        let N_PRIMES = self.n_primes;
        let mut sk : Vec<i32> = vec!();
        for i in 0..N_PRIMES{
            let l_int : i32 = self.l[i].clone().to_str_radix(10).parse().unwrap();
            let m = bound/l_int;
            if m > 0{
                sk.push(rng.gen_range(-m, m) as i32);
            }else{
                sk.push(0);
            }
        }

        let pk = self.class_group_action(Fp::from_int(0), sk.clone());
        (pk, sk)
    }

}

#[cfg(test)]
mod test;