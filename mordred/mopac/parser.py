import re

regexes = [
    ("heat_of_formation", float, re.compile(r"\s*FINAL HEAT OF FORMATION\s*=\s*(-?[0-9.]+)")),
    ("total_energy", float, re.compile(r"\s*TOTAL ENERGY\s*=\s*(-?[0-9.]+)")),
    ("electronic_energy", float, re.compile(r"\s*ELECTRONIC ENERGY\s*=\s*(-?[0-9.]+)")),
    ("core_core_repulsion", float, re.compile(r"\s*CORE-CORE REPULSION\s*=\s*(-?[0-9.]+)")),
    ("ionization_potential", float, re.compile(r"\s*IONIZATION POTENTIAL\s*=\s*(-?[0-9.]+)")),
    ("num_of_filled_levels", int, re.compile(r"\s*NO. OF FILLED LEVELS\s*=\s*(-?[0-9.]+)")),
]
dipole_re = re.compile(r"^ SUM +(?:-?[0-9.]+) +(?:-?[0-9.]+) +(?:-?[0-9.]+) +(-?[0-9.]+)")


class Result(object):
    __slots__ = (
        "heat_of_formation",
        "total_energy",
        "electronic_energy",
        "core_core_repulsion",
        "ionization_potential",
        "num_of_filled_levels",
        "eigenvalues",
        "coordinates",
        "dipole",
    )

    @property
    def occupied_molecular_orbitals(self):
        return self.eigenvalues[:self.num_of_filled_levels]

    @property
    def unoccupied_molecular_orbitals(self):
        return self.eigenvalues[self.num_of_filled_levels:]

    @property
    def homo(self):
        return max(self.occupied_molecular_orbitals)

    @property
    def lumo(self):
        return min(self.unoccupied_molecular_orbitals)

    def to_object(self):
        o = {v: getattr(self, v) for v in self.__slots__ if hasattr(self, v)}
        o["homo"] = self.homo
        o["lumo"] = self.lumo
        return o


def skip_to(s, it):
    for line in it:
        if s in line:
            break


def skip(n, it):
    for _ in range(n):
        next(it)


def parse_output(fp):
    i = 0
    result = Result()

    with fp:
        it = iter(fp)
        for line in it:
            name, conv, regex = regexes[i]
            matched = regex.match(line)
            if matched is not None:
                setattr(result, name, conv(matched.group(1)))
                i += 1
                if i >= len(regexes):
                    break

        skip_to("EIGENVALUES", it)
        skip(1, it)

        evs = []
        for line in it:
            line = line.rstrip()
            if not line:
                break

            for i in range(int(len(line) // 10)):
                evs.append(float(line[10 * i:10 * (i + 1)]))

        result.eigenvalues = evs

        skip_to("DIPOLE           X         Y         Z       TOTAL", it)
        skip(2, it)
        matched = dipole_re.match(next(it))
        if matched is not None:
            result.dipole = float(matched.group(1))

        skip_to("CARTESIAN COORDINATES", it)
        skip(3, it)

        crd = []
        for line in it:
            line = line.rstrip()
            if not line:
                break

            crd.append([float(line[30:40]), float(line[40:50]), float(line[50:60])])

        result.coordinates = crd

    return result


def main():
    import sys
    import json
    json.dump(parse_output(sys.stdin).to_object(), sys.stdout, separators=(",", ":"))


if __name__ == "__main__":
    main()
