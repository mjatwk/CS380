# CS380 Intro to Computer Graphics

>  2023 Spring KAIST CS380 Intro to Computer Graphics (Prof. Minhyuk Sung) 
>
> [Course Website][Website] | [Campuswire][Campuswire] | [Gradescope][Gradescope] 


### TA Help
- For QnA, use [Campuswire][Campuswire] or book an office hour via [youcanbookme][youcanbookme]

### Make & Execute each assignment
```shell
make                      # compile the cod
make clean && make        # recompile the code
./asst{assignment number} # execute the binary file
```

for solution binaries in mac
```shell
chmod +x asst2_{intel,arm}_mac
xattr -d com.apple.quarantine asst2_{intel,arm}_mac
./asst2_{intel,arm}_mac
```

### Incorrect parts of Assignments
- asst2 : the translation direction of the case when the view and manipulated object is identical is incorrect; the camera should elevate while clicking right-click and move upwards, when view and manipulated object is both camera

[Website]: https://mhsung.github.io/kaist-cs380-spring-2023/
[Campuswire]: https://campuswire.com/c/G7A5A8CF5/feed
[Gradescope]: https://www.gradescope.com/courses/515340
[youcanbookme]: https://minhyuksung.youcanbook.me/