# Ported from time in Julia Base in order to print to STDERR instead of STDOUT.

macro timer(ex)
    quote
        local stats = Base.gc_num()
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        local diff = Base.GC_Diff(Base.gc_num(), stats)
        timer_print(elapsedtime, diff.allocd, diff.total_time,
                   Base.gc_alloc_count(diff))
        val
    end
end

function timer_print(elapsedtime, bytes, gctime, allocs)
    @printf(STDERR, "%10.6f seconds", elapsedtime/1e9)
    if bytes != 0 || allocs != 0
        bytes, mb = Base.prettyprint_getunits(bytes, length(Base._mem_units), Int64(1024))
        allocs, ma = Base.prettyprint_getunits(allocs, length(Base._cnt_units), Int64(1000))
        if ma == 1
            @printf(STDERR, " (%d%s allocation%s: ", allocs, Base._cnt_units[ma], allocs==1 ? "" : "s")
        else
            @printf(STDERR, " (%.2f%s allocations: ", allocs, Base._cnt_units[ma])
        end
        if mb == 1
            @printf(STDERR, "%d %s%s", bytes, Base._mem_units[mb], bytes==1 ? "" : "s")
        else
            @printf(STDERR, "%.3f %s", bytes, Base._mem_units[mb])
        end
        if gctime > 0
            @printf(STDERR, ", %.2f%% gc time", 100*gctime/elapsedtime)
        end
        print(STDERR, ")")
    elseif gctime > 0
        @printf(STDERR, ", %.2f%% gc time", 100*gctime/elapsedtime)
    end
    println(STDERR)
end
