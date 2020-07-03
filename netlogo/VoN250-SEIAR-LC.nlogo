turtles-own
[
  infected?           ;; if true, the turtle is infectious
  resistant?          ;; if true, the turtle can't be infected
  virus-check-timer   ;; number of ticks since this turtle's last virus-check
  asymptomatic?
  exposed?
  locked?
  wasLocked?
  class
]

to setup
  clear-all
  reset-ticks
  setup-nodes
  ;; setup-spatially-clustered-network
  ask n-of initial-outbreak-size turtles
    [ become-infected ]
  ask links [ set color white ]
  ;;reset-ticks
end


to setup-nodes

  create-turtles number-of-nodes [ set color red
  set infected? false
  set wasLocked? false
  set class "1a"
  setxy -8 - random 12  12 + random 8
  become-susceptible
  set virus-check-timer random virus-check-frequency]

  ask turtles with [class = "1a"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;;-------
  create-turtles number-of-nodes [ set color yellow
  set infected? false
  set wasLocked? false
  set class "1b"
  setxy 8 + random 12  12 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "1b"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;; 2 -------
  create-turtles number-of-nodes [ set color cyan
  set infected? false
  set wasLocked? false
  set class "2a"
  set shape "box"
  setxy -8 - random 12  4 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "2a"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;;---------
  create-turtles number-of-nodes [ set color magenta
  set infected? false
  set wasLocked? false
  set class "2b"
  set shape "box"
  setxy 8 + random 12  4 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "2b"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;; 3 -------
  create-turtles number-of-nodes [ set color blue
  set infected? false
  set wasLocked? false
  set class "3a"
  set shape "circle"
  setxy -8 - random 12  -4 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "3a"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;;---------
create-turtles number-of-nodes [ set color green
  set infected? false
  set wasLocked? false
  set class "3b"
  set shape "circle"
  setxy 8 + random 12  -4 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "3b"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]

;; 4 --------
create-turtles number-of-nodes [ set color blue
  set infected? false
  set wasLocked? false
  set class "4a"
  set shape "star"
  setxy -8 - random 12  -12 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "4a"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;;---------
create-turtles number-of-nodes [ set color green
  set infected? false
  set wasLocked? false
  set class "4b"
  set shape "star"
  setxy 8 + random 12  -12 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "4b"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]

;; 5 -------
create-turtles number-of-nodes [ set color blue
  set infected? false
  set wasLocked? false
  set class "5a"
  set shape "x"
  setxy -8 - random 12  -20 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "5a"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]
;;---------
create-turtles number-of-nodes [ set color green
  set infected? false
  set wasLocked? false
  set class "5b"
  set shape "x"
  setxy 8 + random 12  -20 + random 8

  become-susceptible
    set virus-check-timer random virus-check-frequency]
  ask turtles with [class = "5b"][
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 99 [
        create-link-with myself
      ]
    ]
  ]

  ;; link fra tutti
  ask turtles [
    ask turtles with [ who > [ who ] of myself ] [
      if random 100 < 4 [
        create-link-with myself
      ]
    ]
  ]

;;  repeat 10
;;  [
;;    layout-spring turtles links 0.3 (world-width / (sqrt 30)) 1
;;  ]
end

to go
  if all? turtles [not infected? and not exposed? and not asymptomatic? and not locked?]
    [ stop ]
  ask turtles
  [

     if locked? [
     set virus-check-timer virus-check-timer + 1
     if virus-check-timer >= v-L-check-freq
       [ set virus-check-timer 0 ]
      ]

    if infected? [
     set virus-check-timer virus-check-timer + 1
     if virus-check-timer >= virus-check-frequency
       [ set virus-check-timer 0 ]
      ]

    if asymptomatic? [
     set virus-check-timer virus-check-timer + 1
     if virus-check-timer >= virus-check-frequency
       [ set virus-check-timer 0 ]
      ]

    if exposed? [
     set virus-check-timer virus-check-timer + 1
     if virus-check-timer >= v-E-check-freq
       [ set virus-check-timer 0 ]
      ]
  ]
  spread-virus
  do-virus-checks
  tick
end

to become-locked  ;; turtle procedure
;;  print "locking!"
  set infected? false
  set asymptomatic? false
  set resistant? false
  set exposed? false
  set locked? true
  set wasLocked? true
  set color green
  set virus-check-timer 1
  print ticks
;;  print " locked!"
end

to become-exposed  ;; turtle procedure
  set infected? false
  set asymptomatic? false
  set resistant? false
  set exposed? true
  set locked? false
  set color yellow
  set virus-check-timer 1
end

to become-infected  ;; turtle procedure
  set infected? true
  set asymptomatic? false
  set resistant? false
  set exposed? false
  set locked? false
  set color red
  set virus-check-timer 1
end

to become-asymptomatic  ;; turtle procedure
  set infected? false
  set asymptomatic? true
  set resistant? false
  set exposed? false
  set locked? false
  set color cyan
  set virus-check-timer 1
end

to become-susceptible  ;; turtle procedure
  set infected? false
  set asymptomatic? false
  set resistant? false
  set exposed? false
  set locked? false
  set color blue
end

to become-resistant  ;; turtle procedure
  set infected? false
  set asymptomatic? false
  set resistant? true
  set exposed? false
  set locked? false
  set color gray
  ask my-links [ set color gray - 2 ]
end

to spread-virus
   ;; logica di locking
    (foreach ["1a" "1b" "2a" "2b" "3a" "3b" "4a" "4b" "5a" "5b"] [
    x -> ask turtles with [infected? or asymptomatic? and class = x ] [
      ifelse count turtles with [ infected? ] + count turtles with [ asymptomatic? ] >= 3
      [
        if not wasLocked? [
          print "going to lock!"
          ask turtles with [ class = x ] [
            become-locked
          ]
        ]
      ]
      [
          ;; un infetto trasmette ai susceptibles con una certa probabilita'
        ask link-neighbors with [not resistant? and not exposed? and not infected? and not asymptomatic? and not locked?][
          if random-float 100 < virus-spread-chance
          [ ;;print "1a else"
            become-exposed ]
        ]

        ask turtles with [asymptomatic?] [
          ask link-neighbors with [not resistant? and not exposed? and not infected? and not asymptomatic? and not locked?][
            if random-float 100 < virus-spread-chance / 2
               [ become-exposed ]
          ]
        ]
      ]
    ]
  ])

  ;; un infetto trasmette ai susceptibles con una certa probabilita'
;;  ask turtles with [infected?] [
;;    ask link-neighbors with [not resistant? and not exposed? and not infected? and not asymptomatic?][
;;      if random-float 100 < virus-spread-chance
;;        [ become-exposed ]
;;    ]
;;  ]

;;  ask turtles with [asymptomatic?] [
;;    ask link-neighbors with [not resistant? and not exposed? and not infected? and not asymptomatic?][
;;      if random-float 100 < virus-spread-chance / 2
;;        [ become-exposed ]
;;    ]
;;  ]
;; di questo non c'e' bisogno perche' un exposed diventa certamente infetto
;; dopo un certo tempo
;;  ask turtles with [exposed?] [
;;    ask link-neighbors with [not resistant? and not infected?][
;;      if random-float 100 < virus-E2I-chance
;;        [ become-infected ]
;;    ]
;;  ]

end

to do-virus-checks
  ;; dopo un certo tempo un exposed diventa o aymptomatic o infected con una certa
  ;; probabilità, entrambi possono infettare
  ask turtles with [locked? and virus-check-timer = 0]
  [
    become-susceptible
    set wasLocked? true
  ]

  ask turtles with [exposed? and virus-check-timer = 0]
  [
    ifelse random 100 < prob-asymptomatic
      [become-asymptomatic]
      [become-infected]
  ]

  ;; dopo un certo tempo un infected diventa resistant
  ask turtles with [infected? or asymptomatic? and virus-check-timer = 0]
  [
      become-resistant
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
324
16
872
565
-1
-1
13.171
1
10
1
1
1
0
0
0
1
-20
20
-20
20
1
1
1
ticks
30.0

SLIDER
8
165
175
198
virus-spread-chance
virus-spread-chance
0.0
10.0
6.8
0.1
1
%
HORIZONTAL

BUTTON
9
120
104
160
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
115
120
210
160
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
8
346
318
585
Network Status
time
% of nodes
0.0
100.0
0.0
100.0
true
true
"" ""
PENS
"susceptible" 1.0 0 -13345367 true "" "plot (count turtles with [not infected? and not resistant? and not exposed? and not asymptomatic? ]) / (count turtles) * 100"
"infected" 1.0 0 -2674135 true "" "plot (count turtles with [infected?]) / (count turtles) * 100"
"resistant" 1.0 0 -7500403 true "" "plot (count turtles with [resistant?]) / (count turtles) * 100"
"exposed" 1.0 0 -1184463 true "" "plot (count turtles with [exposed?]) / (count turtles) * 100"
"asymptomatic" 1.0 0 -11221820 true "" "plot (count turtles with [asymptomatic?]) / (count turtles) * 100"

SLIDER
9
10
174
43
number-of-nodes
number-of-nodes
5
30
25.0
5
1
NIL
HORIZONTAL

SLIDER
8
200
175
233
virus-check-frequency
virus-check-frequency
1
20
6.0
1
1
ticks
HORIZONTAL

SLIDER
9
80
173
113
initial-outbreak-size
initial-outbreak-size
1
number-of-nodes
1.0
1
1
NIL
HORIZONTAL

SLIDER
9
45
173
78
prob-asymptomatic
prob-asymptomatic
0
100
30.0
1
1
NIL
HORIZONTAL

BUTTON
224
121
311
161
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
9
241
140
274
v-E-check-freq
v-E-check-freq
0
10
3.0
1
1
NIL
HORIZONTAL

SLIDER
12
281
138
314
v-L-check-freq
v-L-check-freq
0
20
3.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

Model for the spread of an influenza like disease through a contact network.  
Each node may be in one of the following states:  susceptible, exposed, infected-symptomatic, infected-asymptomatic or resistant.  A modified SEIR model for epidemics.

## HOW IT WORKS

Each time step (tick), each infected node (colored red) attempts to infect all of its neighbors.  Susceptible neighbors (colored blue) will become exposed (yellow) with a probability given by the virus-spread-chance slider after
v-E-check-freq ticks an exposed node becomes asymptomatic with probability (prob-asymptomatic) or infected with probability (1 - prob-asymptomatic).
After virus-check-frequency ticks an infected or asymptomatic node becomes resistant (gray) and cannot be infected.  

If 3 nodes in a class are detected the entire class is locked for v-L-check-freq ticks 
after that delay the class come back to susceptibles.


When a node becomes resistant, the links between it and its neighbors are darkened, since they are no longer possible vectors for spreading the virus.


## References
This simulation in inspired by the model descripted in the following article.

Gemmetto, V., Barrat, A. & Cattuto, C. Mitigation of infectious disease at school: targeted class closure vs school closure. BMC Infect Dis 14, 695 (2014)

and 

* Stonedahl, F. and Wilensky, U. (2008).  NetLogo Virus on a Network model.  http://ccl.northwestern.edu/netlogo/models/VirusonaNetwork.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
