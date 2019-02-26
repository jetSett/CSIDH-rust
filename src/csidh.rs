use std::marker::PhantomData;

use rand::Rng;

use elliptic_curve_algorithms::field::*;

use elliptic_curve_algorithms::finite_fields::FiniteField;

#[macro_use] use elliptic_curve_algorithms;
use elliptic_curve_algorithms::elliptic_curves::*;
use elliptic_curve_algorithms::elliptic_curves::fp_elliptic_curves::*;

pub type Integer = gmp::mpz::Mpz;

type PrimeList = Vec<Integer>;

pub struct CSIDHInstance<K>{
    data : PhantomData<K>,
    pub n_primes : usize,
    pub l: PrimeList,
    pub p: Integer,
}


impl<K : FiniteField<Integer=Integer>> CSIDHInstance<K>{

    pub fn new(n_primes : usize, p: Integer, l: PrimeList) -> CSIDHInstance<K>{
        CSIDHInstance{
            data: PhantomData,
            n_primes,
            l, p
        }
    }

    pub fn check_well_defined(self : &CSIDHInstance<K>){
        let L = &self.l;
        let P = &self.p;
        let N_PRIMES = self.n_primes;

        assert_eq!(N_PRIMES, self.l.len());

        let mut prod = Integer::from(1);
        for i in 0..N_PRIMES{
            prod *= L[i].clone();
        }
        assert_eq!(*P, K::cardinal());
        assert_eq!(*P, Integer::from(4)*prod-Integer::from(1));
    }

    fn is_supersingular(self : &CSIDHInstance<K>, ell : &EllipticCurve<K>) -> bool{
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

    fn order_naive(ell : &EllipticCurve<K>, p : &UnsignedProjPoint<K>) -> Integer{
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

    fn isogeny(ell : &EllipticCurve<K>, point : &UnsignedProjPoint<K>, mut q : UnsignedProjPoint<K>, k : Integer) -> Result<(EllipticCurve<K>, UnsignedProjPoint<K>), ()>{
        assert!(ell.is_montgomery());
        assert!(k>=Integer::from(3));
        assert!(&k%2 != Integer::from(0));

        let p = point.clone();

        let mut i = Integer::from(1);
        
        let mut t = p.clone();  // t = [i]p will iterate over elements of <p>
        let mut t_minus_1 = UnsignedProjPoint::infinite_point();

        let mut pi = UnsignedProjPoint::finite_point(K::from_int(1));
        let mut sigma = K::from_int(0);
        
        let mut projection_x = K::from_int(1);
        let mut projection_z = K::from_int(1);

        while &Integer::from(2)*&i < k{
            if t.x == K::from_int(0){ // point of order 2, we abort
                return Err(());
            }

            pi.x *= t.x.clone();
            pi.z *= t.z.clone();

            sigma += t.x.clone()/t.z.clone() - t.z.clone()/t.x.clone();

            projection_x *= t.x.clone()*q.x.clone() - t.z.clone()*q.z.clone();
            projection_z *= q.x.clone()*t.z.clone() - t.x.clone()*q.z.clone();

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

        sigma *= K::from_int(2);

        pi = pi.normalize();
        Ok((
            EllipticCurve::new_montgomery(pi.x*(ell.a_2.clone() - (K::from_int(3)*sigma))),
            UnsignedProjPoint{
                    x: q.x*projection_x,
                    z: q.z*projection_z
                }
        ))

    }

    pub fn verify_public_key(self : &CSIDHInstance<K>, pk : K) -> bool{
        self.is_supersingular(&EllipticCurve::new_montgomery(pk))
    }

    pub fn naive_class_group_action(self : &CSIDHInstance<K>, pk : K, sk : Vec<i32>) -> K{
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
                let compute_rhs = | x : K | {
                ( (x.clone()*x.clone())*x.clone()  + ((ell.a_2.clone())*x.clone())*x.clone() )+ x 
                };

                let mut x = K::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));
                let mut p_point = UnsignedProjPoint::finite_point(x.clone());

                let mut q_point = ell.scalar_mult_unsigned(p_plus_1.clone()/L[i].clone(), p_point.clone());

                while compute_rhs(x.clone()).legendre_symbol() != s as i8 || q_point == UnsignedProjPoint::infinite_point() {
                    x = K::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));

                    p_point = UnsignedProjPoint::finite_point(x.clone());
                    q_point = ell.scalar_mult_unsigned(p_plus_1.clone()/L[i].clone(), p_point);
                }

                let (ell_, _) = Self::isogeny(&ell, &q_point, q_point.clone(), L[i].clone()).unwrap();
                ell = ell_;
            }

        }
        ell.a_2
    }

    pub fn class_group_action(self : &CSIDHInstance<K>, pk : K, mut sk : Vec<i32>) -> K{
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
            let x = K::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));
            let compute_rhs = | x : K | {
            ( (x.clone()*x.clone())*x.clone()  + ((ell.a_2.clone())*x.clone())*x.clone() )+ x 
            };

            let s = compute_rhs(x.clone()).legendre_symbol() as i32;
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

                let r_point = ell.scalar_mult_unsigned(k.clone()/(li.clone()), q_point.clone());

                if &r_point == &UnsignedProjPoint::infinite_point(){
                    // The point sampled had too low l-index
                    continue;
                }

                let (ell_, q_point_) = match Self::isogeny(&ell, &r_point, q_point.clone(), li.clone()){
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

    pub fn sample_keys(self : &CSIDHInstance<K>, m : i32) -> (K, Vec<i32>){
        let mut rng = rand::thread_rng();
        let N_PRIMES = self.n_primes;
        let mut sk : Vec<i32> = vec!();
        for _i in 0..N_PRIMES{
            sk.push(rng.gen_range(-m, m) as i32);
        }

        let pk = self.class_group_action(K::from_int(0), sk.clone());
        (pk, sk)
    }

}

#[cfg(test)]
mod test;