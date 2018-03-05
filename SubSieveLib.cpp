#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <set>

using namespace std;

#define VERBOSE true

// Activating/DeActivating Tricks tricks

#define PROGRESSIVE_SIEVING true
#define XOR_POPCNT_TRICK true
#define XOR_POPCNT_THRESHOLD 47


typedef float FT;       // Our floating point type
typedef int16_t ZT;     // Our floating integer type
typedef array<uint64_t, 2> CompressedVector;    // Compressed vector type for XOR-POPCNT


// A vector with two representations, and side information
struct Entry {
  vector<ZT> x;           // Vector coordinates in basis B
  vector<FT> y;           // Vector coordinates in gso basis
  vector<FT> yr;          // Vector coordinates in gso basis renormalized (for faster inner product)
  FT l = 0.;              // (squared) length of the vector
  CompressedVector c;     // Compressed vector
  uint64_t uid;           // Unique identifier for collision detection
};

// A compressed vector, and pointer to the full vector
struct CompressedEntry {
  size_t i;               // Index of the non compressed entry
  CompressedVector c;     // Store the compressed vector here to avoid an indirection
  FT l;                   // Sorting Value
};


// Overloaded operations
bool compare_CompressedEntry(CompressedEntry const& lhs, CompressedEntry const& rhs) 
{
    return lhs.l < rhs.l;
}


// Main class
class Siever
{
public:
    int n;                           // current dimension (correspiond to n-d in the paper, until call to the lift)
    vector<vector<FT>> mu;           // gso coefficients, triangular matrix
    vector<FT> r;                    // gso squared norms (normalized by gh)

    vector<Entry> db;                // database
    vector<CompressedEntry> cdb;     // compressed version, faster access and sortable
    unordered_set<uint64_t> db_uid;

    unsigned compress_pos[128][4];   // Indices to chose for compression
    vector<uint64_t> uid_coeffs;
    unordered_set<size_t> updated;
    long unsigned int collisions = 0;
    vector<Entry> projectors;
    int remplace_index;


    // initialize the siever, loading a sub-GSO of dimension nn-d
    void initialize(int nn, int d, double* mu_e, double gh_e)
    {
        n = nn - d;
        fill_compress_pos();
        fill_uid_coeffs();

        mu.resize(n);
        r.resize(n);
        if (VERBOSE) cerr << "Siever: dimension set to " << n << endl;
        for (int i = 0; i < n; ++i)
        {
            mu[i].resize(n);
            for (int j = 0; j < n; ++j)
            {      
                mu[i][j] = mu_e[(i+d)*nn + (j+d)];
                //cerr << mu[i][j] << " ";
            }
            //cerr << endl;
        }
        for (int i = 0; i < n; ++i)
        {
            r[i] = mu[i][i] / gh_e;
            mu[i][i] = 1;
        }
        if (VERBOSE) cerr << "gso loaded  " << n << endl;
        if (VERBOSE) cerr << "gh  " << gh_e << endl;
    }


    // Increase the GSO to the left from n-d to n, and lift the database using babai
    void lift(int nn, int d, double* mu_e, double gh_e)
    {
        initialize(nn, 0, mu_e, gh_e);

        uid_coeffs.clear();
        for (CompressedEntry& ce : cdb)
        {
            for (int i = 0; i < d; ++i)
            {                
                db[ce.i].x.insert(db[ce.i].x.begin(), 0);
            }
            db[ce.i].y.resize(n);
            db[ce.i].yr.resize(n);

            to_gs(db[ce.i].y, db[ce.i].x);
            babai(db[ce.i].y, db[ce.i].x, 0, d);
            to_gs(db[ce.i].y, db[ce.i].x);
            normalize(db[ce.i].yr, db[ce.i].y);
            db[ce.i].l = length(db[ce.i].y);

            db[ce.i].uid = compute_uid(db[ce.i].x);
            db[ce.i].c = compress(db[ce.i].yr);
            ce.l = db[ce.i].l;
            ce.c = db[ce.i].c;
            db_uid.insert(db[ce.i].uid);
            updated.insert(&ce - &cdb[0]);
        }        

        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);
    }

    // Insert a projector for the attempted HKZ reduction
    void insert_projector(long* x)
    {
        uid_coeffs.clear();
        Entry e;
        e.x.resize(n);
        e.y.resize(n);
        for (int i = 0; i < n; ++i)
        {
            e.x[i] = x[i];
        }
        to_gs(e.y, e.x);
        e.l = length(e.y);
        projectors.push_back(e);


        for (CompressedEntry& ce : cdb)
        {

            to_gs(db[ce.i].y, db[ce.i].x);
            db[ce.i].l = length(db[ce.i].y);
            if (db[ce.i].l < 1e-3)
            {
                for (int j = 0; j < n; ++j)
                {
                    db[ce.i].x[j] = (rand() % 3) - 1;
                }            
                to_gs(db[ce.i].y, db[ce.i].x);
                db[ce.i].l = length(db[ce.i].y);
            }

            normalize(db[ce.i].yr, db[ce.i].y);
            db[ce.i].uid = compute_uid(db[ce.i].x);            
            db[ce.i].c = compress(db[ce.i].yr);
            ce.l = db[ce.i].l;
            ce.c = db[ce.i].c;
            db_uid.insert(db[ce.i].uid);
            updated.insert(&ce - &cdb[0]);
        }
        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);
        if (VERBOSE)
        {
            cerr << "After projection";
            //print_db_stats();
        }
    }

    // Load an external database of size N
    void load_db(int N, long* db_e)
    {
        Entry e;
        e.x.resize(n);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                e.x[j] = db_e[i*n + j];
            }
            insert_in_db(e);
        }
    }

    // Save to an external database (select the N shortest ones)
    void save_db(int N, long* db_e)
    {
        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);
        if (N > cdb.size())
        {
            cerr << " Asking for too much in save_db" << endl;
            exit(1);
        }

        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                db_e[i*n + j] = db[cdb[i].i].x[j];
            }           
        }
    }

    // Convert a integer vector in base B to a real vector in base B*
    // and then apply projections from projectors if any
    inline void to_gs(vector<FT> &y,const vector<ZT> &x)
    {
        FT l = 0;
        fill(y.begin(), y.end(), 0.);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                y[j] += x[i] * mu[i][j];
            }
        }
        for (Entry& e : projectors)
        {
            FT inner = 0;
            for (int i = 0; i < n; ++i)
            {
                inner += y[i] * e.y[i] * r[i];
            }
            inner /= e.l;
            for (int i = 0; i < n; ++i)
            {
                y[i] -= inner * e.y[i];
            }            
        }
    }

    // Apply babai algorithm from indices n2 downto n1
    inline void babai(vector<FT> &y, vector<ZT> &x, int n1, int n2)
    {
        int c;
        for (int i = n2-1; i >= n1; --i)
        {
            c = - floor(y[i] + 0.5);
            x[i] += c;
            for (int j = 0; j <= i; ++j)
            {
                y[j] += c * mu[i][j];
            }
        }
    }

    // Compute the renormalized representation for later faster inner products
    inline void normalize(vector<FT> &yr,const vector<FT> &y)
    {
        for (int i = 0; i < n; ++i)
        {
            yr[i] = y[i] * sqrt(r[i]);
        }
    }

    // compute the length of a vector
    inline FT length(vector<FT> &y)
    {
        FT l = 0.;
        for (int i = 0; i < n; ++i)
        {
            l += y[i] * y[i] * r[i];
        }
        return l;
    }

    // choose coefficients for the unique identifier function
    void fill_uid_coeffs()
    {
        uid_coeffs.resize(n);
        for (unsigned i = 0; i < n; i++)
        {
            uid_coeffs[i] = (((uint64_t) rand()) << 32) + ((uint64_t) rand());
        }
        return;
    }

    // compute unique identifier for vectors (invariant by symmetry)
    long unsigned int inline compute_uid(vector<ZT> &x)
    {
        uint64_t uid = 0;
        long int sign = 0;
        for (int i = 0; i < n; ++i)
        {
            if (!sign)
            {
                sign += x[i] > 0;
                sign -= x[i] < 0;
            }
            uid += sign * x[i] * uid_coeffs[i];            
        }
        return uid;
    }

    // choose the vectors sparse vectors r_i for the compressed representation
    void fill_compress_pos()
    {
        for (unsigned i = 0; i < 128; i++)
        {
            compress_pos[i][0] = i % n;
            compress_pos[i][1] = rand() % n;
            compress_pos[i][2] = rand() % n;
            compress_pos[i][3] = rand() % n;
        }
    }


    // Compute the compressed representation of an entry
    CompressedVector inline compress(vector<FT> &v)
    {
        CompressedVector c = {0, 0};
        uint64_t clo, chi;
        float a[2];
        uint32_t *alo = (uint32_t *)(&a[0]);
        uint32_t *ahi = (uint32_t *)(&a[1]);

        clo = 0;
        chi = 0;
        for (int i = 0; i < 64; i++)
        {
            a[0]   = v[compress_pos[2*i][0]];
            a[0]  -= v[compress_pos[2*i][1]];
            a[0]  += v[compress_pos[2*i][2]];
            a[0]  -= v[compress_pos[2*i][3]];

            a[1]   = v[compress_pos[2*i + 1][0]];
            a[1]  -= v[compress_pos[2*i + 1][1]];
            a[1]  += v[compress_pos[2*i + 1][2]];
            a[1]  -= v[compress_pos[2*i + 1][3]];

            clo |= (uint64_t)((alo[0] & 0x80000000) >> 31) << i;
            chi |= (uint64_t)((ahi[0] & 0x80000000) >> 31) << i;
        } 
        c[0] = clo;
        c[1] = chi;
        return c;
    }

    // insert an entry in the database
    bool insert_in_db(Entry &e)
    {
        e.y.resize(n);
        e.yr.resize(n);
        e.uid = compute_uid(e.x);
        unordered_set<uint64_t>::iterator got;
        got = db_uid.find(e.uid);
        if (got != db_uid.end()) 
        {
            collisions ++;
            return false;
        }
        to_gs(e.y, e.x);
        normalize(e.yr,e.y);
        e.c = compress(e.yr);
        e.l = length(e.y);
        if (e.l < 1e-3) { return false;}
        db_uid.insert(e.uid);
        db.push_back(e);
        CompressedEntry ce;
        ce.i = db.size()-1;
        ce.l = e.l;
        ce.c = e.c;
        cdb.push_back(ce);
        updated.insert(cdb.size()-1);
        return true;
    }

    // remplace an entry in the database
    bool remplace_in_db(CompressedEntry *ce, Entry &e)
    {
        e.uid = compute_uid(e.x);
        unordered_set<uint64_t>::const_iterator got;
        got = db_uid.find(e.uid);
        if (got != db_uid.end()) 
        {
            collisions++;
            //cerr << "C";
            return false;
        }
        to_gs(e.y, e.x);
        normalize(e.yr,e.y);
        e.c = compress(e.yr);
        e.l = length(e.y);
        if (e.l < 1e-3) 
        { 
            //cerr << "L"; 
            return false;
        }
        if (1.001 * e.l > ce->l) 
        { 
            //cerr << "Z";
            return false;
        }
        db_uid.erase(db[ce->i].uid);
        db_uid.insert(e.uid);
        db[ce->i] = e;
        ce->l = e.l;
        ce->c = e.c;
        updated.insert(ce - &cdb[0]);
        //cerr << "P";
        return true;
    }

    // Sample a fresh vector, usiong only the partial dimension n1
    bool sample_and_insert_in_db(int n1)
    {
        Entry e;
        e.x.resize(n);
        e.y.resize(n);
        bool res = false;
        while(!res){
            fill(e.x.begin(), e.x.end(), 0);
            int n2 =  n1 / 2;
            for (int j = n2; j < n1; ++j){
                e.x[j] = (rand() % 5) - 2;
            }

            to_gs(e.y, e.x);
            babai(e.y, e.x, 0, n2);

            res = insert_in_db(e);

        }
        return res;
    }

    inline bool reduce(CompressedEntry *ce1, CompressedEntry *ce2)
    {
        FT inner = 0.;
        FT* yr1 = & (db[ce1->i].yr.front());
        FT* yr2 = & (db[ce2->i].yr.front());

        for (int k = 0; k < n; ++k)
        {
            inner += yr1[k] * yr2[k];
        }
        FT new_l = 1.001 * (ce1->l + ce2->l - 2 * abs(inner));
        if ((new_l > ce1->l) && (new_l > ce2->l))
        {
            return false;
        }
        if ((ce1->l < 1e-3) || (ce2->l < 1e-3))
        {
            return false;
        }
        int sign = inner > 0 ? 1 : -1;
        Entry new_e;
        new_e.x.resize(n);
        new_e.y.resize(n);
        new_e.yr.resize(n);
        for (int j = 0; j < n; ++j)
        {
            new_e.x[j] = db[ce1->i].x[j] - sign * db[ce2->i].x[j];
        }
        if (ce1->l < ce2->l)
            return remplace_in_db(ce2, new_e);
        else
            return remplace_in_db(ce1, new_e);
    }

/*    // Attempt reduction
    inline bool reduce(CompressedEntry *ce1, CompressedEntry *ce2)
    {
        FT inner = 0.;
        FT* yr1 = & (db[ce1->i].yr.front());
        FT* yr2 = & (db[ce2->i].yr.front());

        for (int k = 0; k < n; ++k)
        {
            inner += yr1[k] * yr2[k];
        }
        FT new_l = 1.001 * (ce1->l + ce2->l - 2 * abs(inner));
        if (new_l > cdb[remplace_index].l)
        {
            return false;
        }
        if ((ce1->l < 1e-3) || (ce2->l < 1e-3))
        {
            return false;
        }
        int sign = inner > 0 ? 1 : -1;
        Entry new_e;
        new_e.x.resize(n);
        new_e.y.resize(n);
        new_e.yr.resize(n);
        for (int j = 0; j < n; ++j)
        {
            new_e.x[j] = db[ce1->i].x[j] - sign * db[ce2->i].x[j];
        }
        if (remplace_in_db(&cdb[remplace_index], new_e))
        {
            remplace_index --;
            return true;
        }
        return false;
    }*/

    // Sieve the current database
    bool sieve_loop()
    {
        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);
        size_t S = cdb.size();
        remplace_index = S - 1;
        //print_db_stats();
        for (int i = 0; i < S; ++i)
        {
            CompressedEntry *pce1 = &cdb[i];
            CompressedEntry* fast_cdb = &cdb[0];
            CompressedVector cv = pce1->c;
            for (int j = S-1; j > i; --j)
            {
                int w = __builtin_popcountl(cv[0] ^ fast_cdb[j].c[0]);
                w    += __builtin_popcountl(cv[1] ^ fast_cdb[j].c[1]);
                if (w< XOR_POPCNT_THRESHOLD || w > (128 - XOR_POPCNT_THRESHOLD))
                {
                    reduce(pce1, &cdb[j]);
                    if (remplace_index < S/2) 
                    {
                        // cerr << "saturated sieve" << endl;
                        return true;
                    }
                }
            }            
        }
        // cerr << "unsaturated sieve" << endl;
        return false;     
    }

    // Decide wether sieving is over
    bool sieve_over(int n1)
    {
        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);
        long count = 0;
        double gh = (1.*n1) / n;
        for (int i = n1; i < n; ++i)
        {
            gh *= pow(r[i], -1. / n1);
        }

        for (CompressedEntry& ce : cdb)
        {
            // cerr << ce.l / gh << "\t";
            if (ce.l > 1.3333 * gh) break;
            count++;
        }
        // cerr << endl;
        // cerr << count << " >? " << .25 * pow(1.3333, n1/2.) << endl;
        return 2. * count > .5 * pow(1.3333, n1/2.);
    }


    // Sieve and fill until saturation of the 4/3*gh ball
    int sieve()
    {
        int n1 = n - 20;
        long S;
        while(true)
        {
            S = 4 * pow(1.333, n1/2.);
            while(cdb.size() < S)
            {
                sample_and_insert_in_db(n1);
            }
            sieve_loop();
            n1 += sieve_over(n1);
            if (n1>n) break;
        }
        print_db_stats(n);
    }

    // Print minimal statistic on the database (min, average and max squared length)
    void print_db_stats(int n1)
    {
        if (!VERBOSE) return;        
        if (cdb.size()==0)
        {
            cerr << "Empty DB" << endl;
            return;
        }

        double gh = (1.*n1) / n;
        for (int i = n1; i < n; ++i)
        {
            gh *= pow(r[i], -1. / n1);
        }

        sort(cdb.begin(), cdb.end(), &compare_CompressedEntry);

        FT av = 0;
        for (CompressedEntry& ce : cdb)
        {
            av += ce.l;
        }
        cerr << " min: " << cdb.begin()->l / gh;
        cerr << " ave: " << av / cdb.size() / gh;
        cerr << " max: " << cdb.rbegin()->l / gh << endl;
    }

    // Print detailed statistic on the database (saturation rates for the x*gh for various x)
    void print_detailed_db_stats()
    {
        if (!VERBOSE) return;
        int S = db.size();
        uint64_t counts[100] = {0};
        for (int j = 0; j < db.size(); ++j)
        {
            float l = db[j].l;
            for (int i = 99; i >=0; --i)  
            {
                if ((1.+ i/100.) < l) break;
                counts[i]++;
            }
        }
        cerr << "Length Stats :" << endl;
        cerr << "Expected ball radius for 2*2^" << log(S)/log(2) << " vectors :" << pow(2 * S, 2./n) << endl;
        // print_db_stats();
        cerr.precision(3);
        for (int i = 16; i < 80; i+=2)
        {
            float rr = 1+i/100.;
            float ratio = 2* counts[i] / pow(rr, .5*n);
            cerr << rr;
            if (!i) {cerr << ".0";}
            if (! (i%10)) {cerr << "0";}
            cerr << "gh : " << 100*ratio << "\%\t";
            if ((i%16)==14) {cerr << endl;}
        }
        cerr << endl;
    }
}; // End of Siever Class definition




// Making a C library out of the siever class
Siever siever;

extern "C" 
{
    void initialize(int n, int k, double* mu_e, double gh_e)
    {
        siever.initialize(n, k, mu_e, gh_e);
    }

    void sieve(int size_db)
    {
        siever.sieve();
    }

    void lift(int nn, int k, double* mu_e, double gh_e)
    {
        siever.lift(nn, k, mu_e, gh_e);
    }

    void insert_projector(long* y)
    {
        siever.insert_projector(y);
    }

    long db_size()
    {
        return siever.cdb.size();
    }

    long resize_db(long N)
    {
        if (N < siever.cdb.size())
        {
            sort(siever.cdb.begin(), siever.cdb.end(), &compare_CompressedEntry);
            siever.cdb.resize(N);
        }
    }

    void import_db(int N, long* db_e)
    {
        siever.load_db(N, db_e);
    }
    
    void export_db(int N, long* db_e)
    {
        siever.save_db(N, db_e);
    }

}