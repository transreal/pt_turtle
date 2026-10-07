/* turtle_survey.cu -- GPU turtle scan on a rhomb tiling (all tiles CCW).
   Input  (little endian):
     int32 nTiles, nStarts, progLen, maxCycles
     int32 edgeNbr[nTiles*4], edgeEntry[nTiles*4]  (-1 = boundary)
     int32 prog[progLen]   (0=S,1=R,2=L)
     int32 starts[nStarts] (cells; 4 edges each -> 4*nStarts threads)
     float32 cents[nTiles*2]
   Output: JSON summary on stdout; per-thread binary (period,flag,orbit,maxr) to argv[2].
   Compile: nvcc -O3 -o turtle_survey.exe turtle_survey.cu
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>

__global__ void turtleKernel(const int* __restrict__ eN, const int* __restrict__ eE,
                             const int* __restrict__ prog, int progLen, int maxCycles,
                             const int* __restrict__ starts, int nThreads,
                             const float* __restrict__ cents,
                             int* outPeriod, int* outFlag, int* outOrbit, float* outMaxR)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nThreads) return;
    int c0 = starts[tid >> 2], e0 = tid & 3;
    int c = c0, e = e0;
    int orb = c0 * 4 + e0;
    float x0 = cents[2 * c0], y0 = cents[2 * c0 + 1];
    float mr2 = 0.f;
    for (int cyc = 1; cyc <= maxCycles; cyc++) {
        for (int ci = 0; ci < progLen; ci++) {
            int cmd = prog[ci];
            if (cmd == 0) {
                int idx = c * 4 + e;
                int nb = eN[idx];
                if (nb < 0) { outPeriod[tid] = -1; outFlag[tid] = 1; outOrbit[tid] = orb; outMaxR[tid] = sqrtf(mr2); return; }
                e = (eE[idx] + 2) & 3; c = nb;
                float dx = cents[2 * c] - x0, dy = cents[2 * c + 1] - y0;
                float r2 = dx * dx + dy * dy; if (r2 > mr2) mr2 = r2;
            } else if (cmd == 1) { e = (e + 3) & 3; }
            else { e = (e + 1) & 3; }
        }
        int st = c * 4 + e; if (st < orb) orb = st;
        if (c == c0 && e == e0) { outPeriod[tid] = cyc; outFlag[tid] = 0; outOrbit[tid] = orb; outMaxR[tid] = sqrtf(mr2); return; }
    }
    outPeriod[tid] = -1; outFlag[tid] = 2; outOrbit[tid] = orb; outMaxR[tid] = sqrtf(mr2);
}

struct ClassInfo { int period, ncells; long long norbits, nstarts; double maxrad; int repc, repe; double cx, cy; };

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s in.bin out.bin [maxClassesPrint]\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb"); if (!f) { perror("in"); return 1; }
    int nTiles, nStarts, progLen, maxCycles;
    fread(&nTiles, 4, 1, f); fread(&nStarts, 4, 1, f); fread(&progLen, 4, 1, f); fread(&maxCycles, 4, 1, f);
    std::vector<int> eN(nTiles * 4), eE(nTiles * 4), prog(progLen), starts(nStarts);
    std::vector<float> cents(nTiles * 2);
    fread(eN.data(), 4, nTiles * 4, f); fread(eE.data(), 4, nTiles * 4, f);
    fread(prog.data(), 4, progLen, f); fread(starts.data(), 4, nStarts, f);
    fread(cents.data(), 4, nTiles * 2, f); fclose(f);
    int nThreads = nStarts * 4;
    fprintf(stderr, "nTiles=%d nStarts=%d progLen=%d maxCycles=%d threads=%d\n", nTiles, nStarts, progLen, maxCycles, nThreads);

    int *d_eN, *d_eE, *d_p, *d_s, *d_per, *d_fl, *d_orb; float *d_c, *d_mr;
    cudaMalloc(&d_eN, nTiles * 16); cudaMalloc(&d_eE, nTiles * 16); cudaMalloc(&d_p, progLen * 4);
    cudaMalloc(&d_s, nStarts * 4); cudaMalloc(&d_c, nTiles * 8);
    cudaMalloc(&d_per, nThreads * 4); cudaMalloc(&d_fl, nThreads * 4); cudaMalloc(&d_orb, nThreads * 4); cudaMalloc(&d_mr, nThreads * 4);
    cudaMemcpy(d_eN, eN.data(), nTiles * 16, cudaMemcpyHostToDevice);
    cudaMemcpy(d_eE, eE.data(), nTiles * 16, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p, prog.data(), progLen * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(d_s, starts.data(), nStarts * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(d_c, cents.data(), nTiles * 8, cudaMemcpyHostToDevice);
    int bs = 256, gs = (nThreads + bs - 1) / bs;
    turtleKernel<<<gs, bs>>>(d_eN, d_eE, d_p, progLen, maxCycles, d_s, nThreads, d_c, d_per, d_fl, d_orb, d_mr);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) { fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err)); return 2; }
    std::vector<int> per(nThreads), fl(nThreads), orb(nThreads); std::vector<float> mr(nThreads);
    cudaMemcpy(per.data(), d_per, nThreads * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(fl.data(), d_fl, nThreads * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(orb.data(), d_orb, nThreads * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(mr.data(), d_mr, nThreads * 4, cudaMemcpyDeviceToHost);

    long long nClosed = 0, nEsc = 0, nTO = 0; double maxTOr = 0;
    for (int i = 0; i < nThreads; i++) { if (fl[i] == 0) nClosed++; else if (fl[i] == 1) nEsc++; else { nTO++; if (mr[i] > maxTOr) maxTOr = mr[i]; } }

    /* unique orbits among closed threads */
    std::unordered_map<int, int> orbRep; orbRep.reserve(nThreads);
    std::unordered_map<int, long long> orbCnt;
    for (int i = 0; i < nThreads; i++) if (fl[i] == 0) { if (!orbRep.count(orb[i])) orbRep[orb[i]] = i; orbCnt[orb[i]]++; }

    /* per orbit: recompute visited tiles on CPU, radial profile hash */
    std::vector<int> stamp(nTiles, -1); std::vector<int> visited; visited.reserve(1 << 16);
    std::map<std::tuple<int, int, unsigned long long>, ClassInfo> classes;
    int stampId = 0;
    std::vector<int> rp;
    for (auto& kv : orbRep) {
        int i = kv.second; int c0 = starts[i >> 2], e0 = i & 3; int c = c0, e = e0;
        visited.clear(); stamp[c] = stampId; visited.push_back(c);
        int p = per[i];
        for (int cyc = 0; cyc < p; cyc++) for (int ci = 0; ci < progLen; ci++) {
            int cmd = prog[ci];
            if (cmd == 0) { int idx = c * 4 + e; int nb = eN[idx]; e = (eE[idx] + 2) & 3; c = nb; if (stamp[c] != stampId) { stamp[c] = stampId; visited.push_back(c); } }
            else if (cmd == 1) e = (e + 3) & 3; else e = (e + 1) & 3;
        }
        stampId++;
        double cx = 0, cy = 0; for (int t : visited) { cx += cents[2 * t]; cy += cents[2 * t + 1]; }
        cx /= visited.size(); cy /= visited.size();
        rp.clear(); double maxr = 0;
        for (int t : visited) { double dx = cents[2 * t] - cx, dy = cents[2 * t + 1] - cy; double r = sqrt(dx * dx + dy * dy); if (r > maxr) maxr = r; rp.push_back((int)llround(r * 100.0)); }
        std::sort(rp.begin(), rp.end());
        unsigned long long h = 14695981039346656037ULL;
        for (int v : rp) { h ^= (unsigned long long)(unsigned)v; h *= 1099511628211ULL; }
        auto key = std::make_tuple(p, (int)visited.size(), h);
        auto it = classes.find(key);
        if (it == classes.end()) { ClassInfo ci{p, (int)visited.size(), 1, orbCnt[kv.first], maxr, c0, e0, cx, cy}; classes[key] = ci; }
        else { it->second.norbits++; it->second.nstarts += orbCnt[kv.first]; if (maxr > it->second.maxrad) it->second.maxrad = maxr; }
    }

    printf("{\"nTiles\":%d,\"nStarts\":%d,\"threads\":%d,\"closed\":%lld,\"escape\":%lld,\"timeout\":%lld,\"maxTimeoutRadius\":%.3f,\"nOrbits\":%zu,\"nClasses\":%zu,\"classes\":[\n",
           nTiles, nStarts, nThreads, nClosed, nEsc, nTO, maxTOr, orbRep.size(), classes.size());
    bool first = true;
    for (auto& kv : classes) {
        const ClassInfo& c = kv.second;
        printf("%s{\"period\":%d,\"ncells\":%d,\"norbits\":%lld,\"nstarts\":%lld,\"radius\":%.3f,\"rep\":[%d,%d],\"center\":[%.3f,%.3f]}",
               first ? "" : ",\n", c.period, c.ncells, c.norbits, c.nstarts, c.maxrad, c.repc, c.repe, c.cx, c.cy);
        first = false;
    }
    printf("\n]}\n");

    FILE* fo = fopen(argv[2], "wb");
    if (fo) { fwrite(&nThreads, 4, 1, fo); fwrite(per.data(), 4, nThreads, fo); fwrite(fl.data(), 4, nThreads, fo); fwrite(orb.data(), 4, nThreads, fo); fwrite(mr.data(), 4, nThreads, fo); fclose(fo); }
    cudaFree(d_eN); cudaFree(d_eE); cudaFree(d_p); cudaFree(d_s); cudaFree(d_c); cudaFree(d_per); cudaFree(d_fl); cudaFree(d_orb); cudaFree(d_mr);
    return 0;
}
