    /* Dense CUDA evolution; CPU code above retains initialization and all I/O. */
    gpu_phase("initialization");
    DnsimGpu *gpu = dnsim_gpu_init(&M, sizeof M, g.W, g.Ly, g.a, t0);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    long long t = (long long)t0;
    int absorbed = 0;
    while (t < T) {
        while (li < nlog && tl[li] <= t) li++;
        long long next = li < nlog ? tl[li] : T;
        long long event;
        event = gpu_next_periodic(t, every, 0, T);
        if (event < next) next = event;
        if (fs) {
            event = gpu_next_periodic(t, snap_every, snap_from, T);
            if (event < next) next = event;
        }
        if (ckpt_every) {
            event = gpu_next_periodic(t, ckpt_every, 0, T);
            if (event < next) next = event;
        }
        gpu_phase("compute");
        dnsim_gpu_advance(gpu, (u64)next);
        u64 fl;
        dnsim_gpu_download(gpu, g.a, &fl);
        gpu_phase("io");
        t = next;
        int due = (li < nlog && tl[li] == t) || (every && t % every == 0) || t == T;
        if (due) {
            o = measure(&g, ra, ri);
            fprintf(fo, "%lld %.9f %.9f %.9f %.9f\n", t, o.alive / N, o.bonds / (2 * N), fl / N, o.iface / N);
            if (fp) {
                if (pdir) measure_dir(&g, pa, pb, ra, ri);
                fwrite(&t, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp);
            }
            fflush(fo);
            /* Deliberate difference: absorption is inspected only at output times. */
            if (stop_abs && absorbing && o.iface == 0) { absorbed = 1; break; }
        }
        if (fs && t >= snap_from && (t % snap_every == 0 || t == T)) write_snap(fs, &g, t);
        if (ckpt_every && t % ckpt_every == 0 && t < T) {
            fflush(fo);
            if (fp) fflush(fp);
            if (fs) fflush(fs);
            save_atomic(&g, (u64)t, ckpt);
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    dnsim_gpu_finish(gpu);
