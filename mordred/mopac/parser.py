import re


regexes = [
    ("heat_of_formation", float, re.compile(r"\s*FINAL HEAT OF FORMATION\s*=\s*(-?[0-9.]+)")),
    ("total_energy", float, re.compile(r"\s*TOTAL ENERGY\s*=\s*(-?[0-9.]+)")),
    ("electronic_energy", float, re.compile(r"\s*ELECTRONIC ENERGY\s*=\s*(-?[0-9.]+)")),
    ("core_core_repulsion", float, re.compile(r"\s*CORE-CORE REPULSION\s*=\s*(-?[0-9.]+)")),
    ("ionization_potential", float, re.compile(r"\s*IONIZATION POTENTIAL\s*=\s*(-?[0-9.]+)")),
    ("num_of_filled_levels", int, re.compile(r"\s*NO. OF FILLED LEVELS\s*=\s*(-?[0-9.]+)")),
]


def skip_to(s, it):
    for line in it:
        if s in line:
            break


def skip(n, it):
    for _ in range(n):
        next(it)


def parse_output(fp, eigenvalues=True, coordinates=True):
    i = 0
    result = {}

    with fp:
        it = iter(fp)
        for line in it:
            name, conv, regex = regexes[i]
            matched = regex.match(line)
            if matched is not None:
                result[name] = conv(matched.group(1))
                i += 1
                if i >= len(regexes):
                    break

        if eigenvalues:
            skip_to('EIGENVALUES', it)
            skip(1, it)

            evs = []
            for line in it:
                line = line.rstrip()
                if not line:
                    break

                for i in range(int(len(line) // 10)):
                    evs.append(float(line[10 * i:10 * (i + 1)]))

            result['eigenvalues'] = evs

        if coordinates:
            skip_to('CARTESIAN COORDINATES', it)
            skip(3, it)

            crd = []
            for line in it:
                line = line.rstrip()
                if not line:
                    break

                crd.append([float(line[30:40]), float(line[40:50]), float(line[50:60])])

            result['coordinates'] = crd

    return result


re_dipole = re.compile(r"\s*DIPOLE\s*=\s*(-?[0-9.]+)")


def get_dipole_from_arc(fp):
    with fp:
        for line in fp:
            matched = re_dipole.match(line)
            if matched:
                return float(matched.group(1))


def main():
    import sys
    import json
    json.dump(parse_output(sys.stdin), sys.stdout, separators=(',', ':'))


if __name__ == "__main__":
    main()
