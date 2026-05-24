#include "lzz_pX_extra_extra.h"

/** Pretty print for testing purposes */
bool pp_mon(stringstream& ss, zz_p c, long dx, bool p)
{
    if (c == 0) return false;
    if (p) ss << " + ";

    bool pLast = true;
    if (c != 1) ss << c;
    pLast = false;
    
    if (dx != 0) 
    {
        if (pLast) ss << "*";
        pLast = true;
        
        ss << "x";
    }
    if (dx > 1) ss << "^" << dx;

    return true;
}

string pp(const zz_pX& a)
{
    if (a == 0)
        return "0";
    
    stringstream ss;
    bool p = false;
    
    for (long j = 0; j < a.rep.length(); j++)
        p |= pp_mon(ss, a.rep[j], j, p);

    return ss.str();
}

void CompNaive(zz_pX& c, const zz_pX& a, const zz_pX& b)
{
    if (&c == &a || &c == &b)
    {
        zz_pX cp{};
        CompNaive(cp, a, b);
        c = cp;
        return;
    }

    c.zero();
    long d = deg(a);
    
    for (long i = d; i > -1; i--)
    {
        mul(c, c, b);
        add(c, c, a[i]);
    }
}

void CompInterpolate(zz_pX& c, const zz_pX& a, const zz_pX& b)
{
    if (&c == &a || &c == &b)
    {
        zz_pX cp{};
        CompInterpolate(cp, a, b);
        c = cp;
        return;
    }
    
    c.zero();
    long da = deg(a);
    long db = deg(b);
    long dc = da * db;

    zz_p a0 = random_zz_p();
    zz_p b0 = random_zz_p();

    long xtra = 3;

    zz_pX_Multipoint_Geometric evB;
    evB = zz_pX_Multipoint_Geometric(a0, b0, db + xtra);
    Vec<zz_p> val_b;
    Vec<zz_p> val_c;

    evB.evaluate(val_b, b, dc + 1);
    zz_pX_Multipoint_General evC(val_b);
    evC.evaluate(val_c, a);

    evB = zz_pX_Multipoint_Geometric(a0, b0, dc + 1);
    evB.interpolate(c, val_c);
}

void ShiftX(zz_pX& c, const zz_pX& a, zz_p b)
{
    /*
    [x^i]f(x)=a_i * i!
    [x^i]g(x)=a^{d-1-i} / (d-1-i)!
    [x^{d-1-i}]g(x)=a^i / i!

    and you get [x^i]P(x+a) = [x^{d-1+i}](f*g) / i!

    optimized and works very well !!
   */
    long dx = deg(a);
    zz_pX f = a;
    zz_pX g{}; g.rep.SetLength(a.rep.length());

    zz_p b_to_i(1);
    for (long i = 0; i <= dx; i++)
    {
        f.rep[i] *= fac(i);
        g.rep[dx - i] = b_to_i * inv_fac(i);
        b_to_i *= b;
    }    

    zz_pX h{};
    middle_product(h, f, g, dx, dx);

    c.SetLength(dx + 1);
    for (long i = 0; i <= dx; i++)
        c.rep[i] = h.rep[i] * inv_fac(i); 
    
    c.normalize();
}

void ShiftXInterpolate(zz_pX& c, const zz_pX& a, zz_p b)
{
    if (&c == &a)
    {
        zz_pX cp{};
        ShiftX(cp, a, b);
        c = cp;
        return;
    }
    
    c.zero();
    long da = deg(a);

    zz_p a0 = random_zz_p();
    zz_p b0 = random_zz_p();

    zz_p a02 = a0*a0;
    zz_p pr(1);

    Vec<zz_p> seq;
    seq.SetLength(da + 1);
    for (long i = 0; i < da + 1; i++)
    {
        seq[i] = b0*pr + b;
        pr *= a02;
    }
    zz_pX_Multipoint_General evA(seq);

    Vec<zz_p> val_c;
    evA.evaluate(val_c, a);
    
    zz_pX_Multipoint_Geometric evC(a0, b0, da + 1);
    evC.interpolate(c, val_c);
}
 
zz_pXComp::zz_pXComp(const zz_pX& b, long max_deg_a)
 : a0(random_zz_p()), b0(random_zz_p()), db(deg(b)), evg(a0, b0, db + 1) {
    evg.evaluate(val_b, b, max_deg_a * db + 1);
 }

void eval_partial(vec_zz_p& b, const zz_pX& f, const vec_zz_p& a, long m)
{
    b.SetLength(m);
    long i;
    for (i = 0; i < m; i++) 
        eval(b[i], f, a[i]);
}

void CompInterpolate(zz_pX& c, const zz_pX& a, const zz_pXComp& b)
{
    if (&c == &a)
    {
        zz_pX cp{};
        CompInterpolate(cp, a, b);
        c = cp;
        return;
    }
    
    c.zero();
    long da = deg(a);
    long dc = da * b.db;

    Vec<zz_p> val_c;
    eval_partial(val_c, a, b.val_b, dc + 1); // TODO: use better preconditioned algo

    zz_pX_Multipoint_Geometric evB(b.a0, b.b0, dc + 1);
    evB.interpolate(c, val_c);
}

vec_zz_p zz_pX_roots(const zz_pX& P){
    zz_pX Q_;
    SetCoeff(Q_, 1, to_zz_p(1));

    zz_pX Q;
    PowerMod(Q, Q_, zz_p::modulus(), P);
    Q -= Q_;

    zz_pX RR;
    GCD(RR, P, Q);

    vec_zz_p roots;
    FindRoots(roots, RR);

    return roots;
}