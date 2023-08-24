file_name = "christmas_bear_original"
in_file_path = file_name + ".off"
out_file_path = file_name + ".mesh"
isContainingRGB = True

def main():
    in_file = open(in_file_path, "r")

    n_tri = 0
    n_quad = 0

    lines = []

    #OFF
    in_file.readline()
    
    #v f e
    line = in_file.readline()
    v, f, e = line.split()

    lines.append("0")
    
    for i in range(int(v)):
        x, y, z = in_file.readline().split()
        lines.append(str(x) + " " + str(y) + " " + str(z) + "\n")

    for i in range(int(f)):
        line = in_file.readline()
        n = int(line[0])
        if n == 3:
            if isContainingRGB:
                n, i1, i2, i3, r, g, b = line.split()
            else:
                n, i1, i2, i3 = line.split()
            lines.append(str(i1) + " " + str(i2) + " " + str(i3) + "\n")
            n_tri = n_tri + 1
        elif n == 4:
            if isContainingRGB:
                n, i1, i2, i3, i4, r, g, b = line.split()
            else:
                n, i1, i2, i3, i4 = line.split()
            lines.append(str(i1) + " " + str(i2) + " " + str(i3) + " " + str(i4) + "\n")
            n_quad = n_quad + 1
        else:
            print("[line " + (i+v+3) + "] parsing error: n value not 3 or 4")

    in_file.close()
    
    lines[0] = str(v) + " " + str(n_tri) + " " + str(n_quad) + "\n"
    out_file = open(out_file_path, "w")
    for line in lines:
        out_file.write(line)
    out_file.close()

    print(out_file_path + " successfully generated!")

main()

