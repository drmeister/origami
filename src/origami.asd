(in-package :asdf-user)

(defsystem "origami"
    :description "Build DNA origami"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net> and Tommaso Banelli and Vincenzo Carnevale"
    :licence "LGPL-3.0"
    :depends-on (:read-csv :cl-json)
    :serial t
    :components
    ((:file "packages")
     (:file "parser")
     (:file "build")
     (:file "graph")
     (:file "graphics")
     (:file "energy")
     ))
