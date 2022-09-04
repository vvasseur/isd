import datetime
from gmpy2 import mpz, bincoef
import re
import subprocess
import sys

NB_BENCHMARK = 100


def proba(N, K, W, L, P, EPS):
    s = mpz(0)
    for p in range(0, P + 2, 2):
        P1 = p // 2
        P2 = p - P1
        N1 = (K + L) // 2
        N2 = K + L - N1

        for i1 in range(P1 + 1):
            for i2 in range(P2 + 1):
                s += bincoef(N1 - EPS, P1 - i1) * bincoef(
                    N2 - EPS, P2 - i2) * bincoef(2 * EPS, i1 + i2) * bincoef(
                        N - K - L, W - p)
    if bincoef(N, W) < 2**(N - K):
        s /= bincoef(N, W)
    else:
        s /= 2**(N - K)
    return s


# Estimated time in microseconds
def estimate_time(command, LW, W, P, L, EPS):
    compile_command = """
    cmake -B build -DDUMER_L={l}L -DDUMER_P={p}L -DDUMER_EPS={eps}L -DDUMER_LW={lw} -DBENCHMARK={nb_benchmark}UL && cmake --build build/ -j
    """.format(l=L, p=P, eps=EPS, lw=LW, nb_benchmark=NB_BENCHMARK)
    comp = subprocess.Popen(compile_command,
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    comp.wait()

    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    first_line = proc.stdout.readline().strip().decode('utf-8')
    matches = re.fullmatch(r"n=(\d+) k=(\d+) w=(\d+)", first_line)
    if not matches:
        sys.exit("Error running command.")
    else:
        N = int(matches.group(1))
        K = int(matches.group(2))
        if W == 0:
            W = int(matches.group(3))
    if W == 0:
        sys.exit("Please choose a target weight: W=...")
    last_line = ""
    for line in proc.stdout:
        last_line = line.strip()
    time = 0
    try:
        time = int(last_line)
    except:
        return 0
    return time / proba(N, K, W, L, P, EPS) / 1e3 / NB_BENCHMARK


if __name__ == '__main__':
    command = ["build/isd"]
    P = 4
    W = 0
    EPS = 10
    L = 12
    LW = 0
    for arg in sys.argv[1:]:
        match = re.fullmatch(r"(W|P|EPS|L)=(\d+)", arg)
        if match and match.group(1) == "W":
            W = int(match.group(2))
        elif match and match.group(1) == "P":
            P = int(match.group(2))
        elif match and match.group(1) == "EPS":
            EPS = int(match.group(2))
        elif match and match.group(1) == "L":
            L = int(match.group(2))
        else:
            command.append(arg)
        match = re.fullmatch(r"LW", arg)
        if match:
            LW = 1

    print("{:>5} {:>5} {:>5} {}".format("P", "L", "EPS", "Est. time"))

    # Hill-climbing optimization
    local_min = False
    times = {}
    while not local_min:
        vs = []
        for dL in [0, -1, 1]:
            if L + dL > 0 and L + dL <= 64:
                for dEPS in [0, -1, 1]:
                    if EPS + dEPS > 0 and EPS + dEPS <= 64:
                        vs.append((P, L + dL, EPS + dEPS))
        for v in vs:
            if v not in times:
                times[v] = estimate_time(command, LW, W, *v)
                str_time = str(datetime.timedelta(microseconds=int(times[v])))
                if times[v] != 0:
                    print("{:5} {:5} {:5} {}".format(*v, str_time))
        if all([times[vs[0]] < times[v] or times[v] == 0 for v in vs[1:]]):
            local_min = True
        else:
            _, vmin = min((times[v], v) for v in vs)
            P, L, EPS = vmin

    print("-DDUMER_P={p}L -DDUMER_L={l}L -DDUMER_EPS={eps}L".format(
        p=P, l=L, eps=EPS))
