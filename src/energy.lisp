(in-package :origami)


;;; --------------------------------------------------
;;;

(defun build-one-molecule (strand trans)
  (let ((molecule (chem:make-molecule)))
    (loop for node in (dna-strand-as-list-of-nodes strand)
       for res = (chem:matter-copy (residue node))
       do (progn
            (chem:apply-transform-to-atoms res trans)
            (chem:add-matter molecule res))
       do (let* ((hbnode (hbond-node node))
                 (hbres (and hbnode (chem:matter-copy (residue hbnode)))))
            (when hbres
              (chem:apply-transform-to-atoms hbres trans)
              (chem:add-matter molecule hbres))))
    molecule))


(defun build-one-aggregate (origami energy-function pos)
  (let ((strands (parts origami)))
    (let ((agg (chem:make-aggregate)))
      (loop for strand in strands
            for id = (order-id strand)
            for trans = (rigid-body-transform id pos)
            do (let ((m (build-one-molecule strand trans)))
                 (chem:add-matter agg m)))
      agg)))


(defun build-pseudo-aggregate (energy-function position)
  (let* ((terms (chem:rigid-body-energy-function-terms energy-function))
         (nonbond-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-nonbond))
                                terms))
         (staple-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-staple))
                               terms))
         (aggregate (chem:make-aggregate)))
    (when nonbond-term
      do (loop for parts in (chem:parts-as-list nonbond-term position)
               for first-part = (first parts)
               for last-part = (car (last parts))
               for first-atom = (chem:make-atom :f :f)
               for last-atom = (chem:make-atom :o :o)
               for residue = (chem:make-residue :dummy)
               for molecule = (chem:make-molecule :dummy)
               do
                  (chem:set-position first-atom (geom:vec (elt first-part 0)
                                                          (elt first-part 1)
                                                          (elt first-part 2)))
                  (chem:set-position last-atom (geom:vec (elt last-part 0)
                                                         (elt last-part 1)
                                                         (elt last-part 2)))
                  (chem:add-matter residue first-atom)
                  (chem:add-matter residue last-atom)
                  (chem:bond-to first-atom last-atom :single-bond)
                  (chem:add-matter molecule residue)
                  (chem:add-matter aggregate molecule)))
    (when staple-term
      do (loop for parts in (chem:parts-as-list staple-term position)
               for first-atom = (chem:make-atom :n :n)
               for last-atom = (chem:make-atom :n :n)
               for residue = (chem:make-residue :dummy)
               for molecule = (chem:make-molecule :dummy)
               do
                  (chem:set-position first-atom (geom:vec (elt parts 0)
                                                          (elt parts 1)
                                                          (elt parts 2)))
                  (chem:set-position last-atom (geom:vec (elt parts 3)
                                                         (elt parts 4)
                                                         (elt parts 5)))
                  (chem:add-matter residue first-atom)
                  (chem:add-matter residue last-atom)
                  (chem:bond-to first-atom last-atom :single-bond)
                  (chem:add-matter molecule residue)
                  (chem:add-matter aggregate molecule)))
    aggregate))

(defun extract-coordinates-from-aggregate (aggregate coordinates)
  (let ((index -1)) ;; start at -1 to use incf
    (cando:do-atoms (a aggregate)
      (let ((pos (chem:get-position a)))
        (setf (aref coordinates (incf index)) (float (geom:vx pos) 1.0s0)
              (aref coordinates (incf index)) (float (geom:vy pos) 1.0s0)
              (aref coordinates (incf index)) (float (geom:vz pos) 1.0s0))))))

(defun energy-term-of-type (energy-function type)
  (let ((terms (chem:rigid-body-energy-function-terms energy-function))
        (term (find-if (lambda (x) (typep x type))
                       terms)))
    term))

  
(defun rigid-body-energy-function-saved-position (energy)
  (let* ((start-pos (make-array (chem:get-nvector-size energy) :element-type 'double-float)))
    (chem:load-coordinates-into-vector energy start-pos)
    start-pos))

(defun normalize-rigid-body-energy-function-position (energy)
  (let* ((start-pos (make-array (chem:get-nvector-size energy) :element-type 'double-float)))
    (chem:load-coordinates-into-vector energy start-pos)
    (chem:rigid-body-energy-function-normalize-position energy start-pos)
    (chem:save-coordinates-from-vector energy start-pos)))

(defun randomize-rigid-body-energy-function-position (energy)
  (flet ((random-pos ()
           (float (- (random 200.0) 100.0) 1.0d0))
         (random-quat ()
           (float (- (random 2.0) 1.0) 1.0d0)))
    (let* ((start-pos (make-array (chem:get-nvector-size energy) :element-type 'double-float)))
      (loop for index from 0 below (chem:get-nvector-size energy) by *seven*
         do (progn
              (setf (elt start-pos (+ index 0)) (random-quat))
              (setf (elt start-pos (+ index 1)) (random-quat))
              (setf (elt start-pos (+ index 2)) (random-quat))
              (setf (elt start-pos (+ index 3)) (random-quat))
              (setf (elt start-pos (+ index 4)) (random-pos))
              (setf (elt start-pos (+ index 5)) (random-pos))
              (setf (elt start-pos (+ index 6)) (random-pos))))
      (chem:rigid-body-energy-function-normalize-position energy start-pos)
      start-pos)))

(defun arrayed-rigid-body-energy-function-position (energy &key (separation 20.0) verbose)
  (let* ((start-pos (make-array (chem:get-nvector-size energy) :element-type 'double-float)))
    (loop for index from 0 below (chem:get-nvector-size energy) by *seven*
          for num from 0
          do (progn
               (setf (elt start-pos (+ index 0)) 0.0d0)
               (setf (elt start-pos (+ index 1)) 0.0d0)
               (setf (elt start-pos (+ index 2)) 0.0d0)
               (setf (elt start-pos (+ index 3)) 1.0d0)
               (setf (elt start-pos (+ index 4)) (float (* num separation) 1.0d0))
               (setf (elt start-pos (+ index 5)) 0.0d0)
               (setf (elt start-pos (+ index 6)) 0.0d0)))
    (chem:rigid-body-energy-function-normalize-position energy start-pos)
    (when verbose (format t "start-pos -> ~a~%" start-pos))
    start-pos))

#+(or)
(defun move-rigid-body-saved-position (energy-function id vec)
  (let ((pos (make-array (chem:get-nvector-size energy-function) :element-type 'double-float))
        (index (* id *seven*)))
    (chem:load-coordinates-into-vector energy-function pos)
    (setf (elt pos (+ index 4)) (float (first vec) 1.0d0))
    (setf (elt pos (+ index 5)) (float (second vec) 1.0d0))
    (setf (elt pos (+ index 6)) (float (third vec) 1.0d0))
    (chem:save-coordinates-from-vector energy-function pos)))

(defvar *seven* 7)
(defparameter *move-size* 50.0)
(defparameter *rotate-size* 10.0)

(defun jostle-origami (pos &optional (step-size 1.0))
  "Randomly adjust the position and quaternion for the entire origami"
  (let ((half-step (* step-size 0.5)))
    (loop for index from 0 below (length pos) by *seven*
          do (incf (elt pos (+ index 0)) (* *rotate-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 1)) (* *rotate-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 2)) (* *rotate-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 3)) (* *rotate-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 4)) (* *move-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 5)) (* *move-size* (- (random step-size) half-step)))
             (incf (elt pos (+ index 6)) (* *move-size* (- (random step-size) half-step))))))


               

#+(or)
(defun move-rigid-body (energy-function id vec)
  (let ((pos (make-array (chem:get-nvector-size energy-function) :element-type 'double-float))
        (index (* id *seven*)))
    (chem:load-coordinates-into-vector energy-function pos)
    (setf (elt pos (+ index 4)) (float (first vec) 1.0d0))
    (setf (elt pos (+ index 5)) (float (second vec) 1.0d0))
    (setf (elt pos (+ index 6)) (float (third vec) 1.0d0))
    (chem:save-coordinates-from-vector energy-function pos)))

    
(defun build-rigid-body-energy-function (strands)
  (let* ((energy (chem:make-rigid-body-energy-function (length strands))))
    (randomize-rigid-body-energy-function-position energy)
    energy))


(defun check-quaternion (energy-function id)
  (multiple-value-bind (a b c d x y z)
      (chem:rigid-body-energy-function-get-position energy-function id)
    (let ((q (+ (* a a) (* b b) (* c c) (* d d)))
          (m (geom:make-matrix nil)))
      (geom:set-from-quaternion m a b c d x y z)
      (values m q (list a b c d)))))


(defun matrix-from-quaternion (a b c d x y z)
  (let ((m (geom:make-matrix nil)))
    (geom:set-from-quaternion m a b c d x y z)
    m))


(defun add-p-o-staple-term (ef id1 pos1 id2 pos2 &key verbose)
  (when verbose (format t "Adding staple term ef ~a id1 ~a id2 ~a ~%" ef id1 id2 pos1 pos2))
  (chem:energy-rigid-body-staple-add-term ef
                                          230.0 ; 230.0 ; OS-P from parm99.dat
                                          1.61 ; OS-p from parm99.dat
                                          (* *seven* id1)
                                          pos1
                                          (* *seven* id2)
                                          pos2))

(defun build-energy-rigid-body-staple (strands &key verbose)
  (let ((ef (chem:make-energy-rigid-body-staple)))
    (loop for strand in strands
       for p5 = (p5-end strand)
       for p3 = (p3-end strand)
       for p5hb = (hbond-node p5)
       for p3hb = (hbond-node p3)
       for my-order-id = (order-id strand)
       do (let ((p5-partner (backward-node p5)))
            (when (and p5-partner (<= my-order-id (order-id (strand p5-partner))))
              (add-p-o-staple-term ef
                                   my-order-id
                                   (chem:get-position (chem:atom-with-name (residue p5) :P))
                                   (order-id (strand p5-partner))
                                   (chem:get-position (chem:atom-with-name (residue p5-partner) :O3*))
                                   :verbose verbose)))
       do (let ((p3-partner (forward-node p3)))
            (when (and p3-partner (<= my-order-id (order-id (strand p3-partner))))
              (add-p-o-staple-term ef
                                   my-order-id
                                   (chem:get-position (chem:atom-with-name (residue p3) :O3*))
                                   (order-id (strand p3-partner))
                                   (chem:get-position (chem:atom-with-name (residue p3-partner) :P))
                                   :verbose verbose)))
       do (when p5hb
            (let ((p5hb-partner (forward-node p5hb)))
              (when (and p5hb-partner (<= my-order-id (order-id (strand p5hb-partner))))
                (add-p-o-staple-term ef
                                     my-order-id
                                     (chem:get-position (chem:atom-with-name (residue p5hb) :O3*))
                                     (order-id (strand p5hb-partner))
                                     (chem:get-position (chem:atom-with-name (residue p5hb-partner) :P))
                                     :verbose verbose))))
       do (when p3hb
            (let ((p3hb-partner (backward-node p3hb)))
              (when (and p3hb-partner (<= my-order-id (order-id (strand p3hb-partner))))
                (add-p-o-staple-term ef
                                     my-order-id
                                     (chem:get-position (chem:atom-with-name (residue p3hb) :P))
                                     (order-id (strand p3hb-partner))
                                     (chem:get-position (chem:atom-with-name (residue p3hb-partner) :O3*))
                                     :verbose verbose)))))
    ef))


(defun number-of-nonbond-atoms (strand)
  (loop for node in (dna-strand-as-list-of-nodes strand)
     sum (+ (chem:number-of-atoms (residue node))
            (if (hbond-node node)
                (chem:number-of-atoms (residue (hbond-node node)))
                0))))

(defun number-of-nonbond-atoms (strand)
  (let ((total (loop for node in (dna-strand-as-list-of-nodes strand)
                  sum (+ (chem:number-of-atoms (residue node))
                         (if (hbond-node node)
                             (chem:number-of-atoms (residue (hbond-node node)))
                             0)))))
    (format t "Number of nonbond atoms ~a~%" total)
    total))


(defun create-nonbond-terms (nonbond-component strand index)
  (loop for node in (dna-strand-as-list-of-nodes strand)
     for residue = (residue node)
     do (chem:map-atoms
         nil
         (lambda (atom)
           (chem:energy-rigid-body-nonbond-set-term
            nonbond-component
            (prog1 index
              (incf index))
            atom              ; treat everything like a carbon for now
            1.908             ; CT from parm99.dat
            0.1094            ; CT from parm99.dat
            0.0               ; no charge
            (chem:get-position atom)))
                        residue))
  (when (typep strand 'double-stranded-dna)
    (loop for node in (dna-strand-as-list-of-nodes strand)
       for hbnode = (hbond-node node)
       for residue = (residue hbnode)
       do (chem:map-atoms
           nil
           (lambda (atom)
             (chem:energy-rigid-body-nonbond-set-term
              nonbond-component
              (prog1 index
                (incf index))
              atom            ; treat everything like a carbon for now
              1.908           ; CT from parm99.dat
              0.1094          ; CT from parm99.dat
              0.0             ; no charge
              (chem:get-position atom)))
                          residue)))
  index)
       
(defun build-energy-rigid-body-nonbond (strands)
  (format t "strands -> ~a  length -> ~a~%" strands (length strands))
  (let* ((strand-vec (let ((vec (make-array (length strands) :element-type 't))) ; store strands in vector by id
                       (loop for strand in strands ; 
                          do (setf (elt vec (order-id strand)) strand))
                       vec)))
    (format t "(length strand-vec): ~a~%" (length strand-vec))
    (multiple-value-bind (end-atom-vec total-atoms)
                                        ; vector to store the last index of each rigid body
        (let ((vec (make-array (length strands) :element-type 'ext:byte32))
              (total-count 0))
          (loop for id from 0 below (length strand-vec)
             for strand = (elt strand-vec id)
             for strand-atoms-count = (number-of-nonbond-atoms strand)
             for running-count = strand-atoms-count then (+ running-count strand-atoms-count)
             do (setf total-count running-count)
             do (format t "Looking at strand id ~a~%" id)
             do (setf (elt vec id) total-count))
          (values vec total-count))
      (format t "total-atoms: ~a~%" total-atoms)
      (let ((nonbond-component (chem:make-energy-rigid-body-nonbond end-atom-vec))
            (start-index 0))
        (loop for id from 0 below (length strand-vec)
           for strand = (elt strand-vec id)
             do (format t "start-index: ~a ~%" start-index)
           do (setf start-index (create-nonbond-terms nonbond-component strand start-index)))
        nonbond-component))))


(defvar *max-zstep* 10.0)

(defun nonbond-node-z-offset-and-length (strand &key verbose)
  "Calculate the number of nonbond spheres that will represent this dna helix.
Note: This is not the number of bases - the sphere separation *max-zstep* is
larger than the rise between DNA base pairs."
  (let* ((nodes (dna-strand-as-list-of-nodes strand))
         (num-nodes (length nodes)))
    (when (< num-nodes 2)
      (error "A segment of DNA has less than 2 nodes - we can't model this"))
    (let* ((total-height (* (1- num-nodes) +bdna-rise+)) ; Subtract 1 to get the number of steps between nodes
           (z-offset (* total-height -0.5))
           (znum (ceiling (/ total-height *max-zstep*)))
           (zstep (if (/= znum 0)
                      (/ total-height znum)
                      (error "About to divide by zero (znum) z-offset: ~a  total-height: ~a  num-nodes: ~a~%" z-offset total-height num-nodes))))
      (when verbose
        (format t "num-nodes: ~a  total-height: ~a  z-offset: ~a  znum: ~a  zstep: ~a~%"
                num-nodes total-height z-offset znum zstep ))
      (values z-offset total-height znum zstep))))
  
#| Mathematica code to plot enonbond function
      Clear[rab];                       ; ;
      a = { ax, ay, az};                ; ;
      b = {bx, by, bz};                 ; ;
      rab[a_, b_] := Sqrt[(a - b).(a - b)]
      enonbond[a_, b_, r0_, edep_] := 
      edep*((r0/rab[a, b])^12 - 2*(r0/rab[a, b])^6)
      enonbond[{0, 0, 0}, {10, 0, 0}, 1.0, 0.1]
      Plot[enonbond[{0, 0, 0}, {x, 0, 0}, 1.0, 0.01], {x, 0.8, 5}, 
      PlotRange -> Full]
      |#
(defun create-fast-nonbond-terms (nonbond-component strand index offset-start offset-end &key verbose)
  (let* ((nodes (dna-strand-as-list-of-nodes strand))
         (num-nodes (length nodes)))
    (multiple-value-bind (start-height total-height znum zstep)
        (nonbond-node-z-offset-and-length strand)
      (when (< znum 1)
        (error "There is a section of DNA that has ~a nonbond nodes - we can't currently model anything with less than 1" znum))
      (when verbose
        (format t "          Nodes: ~a~%" nodes)
        (format t "      Num nodes: ~a~%" num-nodes)
        (format t "   total-height: ~a~%" total-height)
        (format t "   start-height: ~a~%" start-height)
        (format t "           znum: ~a~%" znum)
        (format t "          zstep: ~a ~%" zstep))
      (loop for node-index from offset-start to (+ znum offset-end)
            for z-height = (+ start-height (* node-index zstep))
            for pos = (geom:vec 0.0 0.0 z-height)
            do (when verbose
                 (format t "node ~a index(~a) z-height ~a~%" node-index index z-height)
                 (format t "Placing at ~a~%" pos))
            do (chem:energy-rigid-body-nonbond-set-term
                nonbond-component
                (prog1 index (incf index))
                (chem:make-atom :bogus :c) ; no atom
                12.0 ; 15 Angstrom radius for a B-form DNA double helix
                10.0 ; Very little attraction
                0.0  ; -2 charge
                pos))
      index)))


(defun build-fast-energy-rigid-body-nonbond (strands &key verbose)
  "Build an energy-rigid-body-nonbond for the strands."
  (when verbose (format t "nonbond term strands -> ~a  length -> ~a~%" strands (length strands)))
  (let ((strand-vec (let ((vec (make-array (length strands) :element-type 't))) ; store strands in vector by id
                      (loop for strand in strands ; 
                            do (setf (elt vec (order-id strand)) strand))
                      vec))
        (offset-start 1) ; start the bottom of the helix one "up"
        (offset-end -1)) ; end the top of the helix one from the top
    (when verbose (format t "(length strand-vec): ~a~%" (length strand-vec)))
    (multiple-value-bind (end-atom-vec total-atoms)
                                        ; vector to store the last index of each rigid body
        (let ((vec (make-array (length strands) :element-type 'ext:byte32))
              (total-count 0))
          (loop for id from 0 below (length strand-vec)
                for strand = (elt strand-vec id)
                for strand-nonbonds-count = (multiple-value-bind (start-height total-height znum zstep)
                                                (nonbond-node-z-offset-and-length strand :verbose verbose)
                                              (- (+ 1 (+ znum offset-end)) offset-start))
                for running-count = strand-nonbonds-count then (+ running-count strand-nonbonds-count)
                do (setf total-count running-count)
                   (setf (elt vec id) running-count)
                   (when verbose (format t "Looking at strand id ~a  nodes: ~a  running-count: ~a~%" id strand-nonbonds-count running-count)))
          (values vec total-count))
      (when verbose
        (format t "total-atoms: ~a~%" total-atoms)
        (format t "end-atom-vec: ~a~%" end-atom-vec))
      (let ((nonbond-component (chem:make-energy-rigid-body-nonbond end-atom-vec))
            (end-index 0))
        (loop for id from 0 below (length strand-vec)
              for strand = (elt strand-vec id)
              do (setf end-index (create-fast-nonbond-terms nonbond-component strand end-index offset-start offset-end :verbose verbose))
              do (when verbose (format t "end-index: ~a ~%" end-index)))
        nonbond-component))))

(defun build-energy-function (origami &key verbose)
  (let ((parts (parts origami)))
    (let ((staple (build-energy-rigid-body-staple parts))
          (nonbonded (build-fast-energy-rigid-body-nonbond parts :verbose verbose))
          (energy (build-rigid-body-energy-function parts)))
      (chem:rigid-body-energy-function-add-term energy staple)
      (chem:rigid-body-energy-function-add-term energy nonbonded)
      energy)))

(defun staple-term (energy-fn)
  (let* ((terms (chem:rigid-body-energy-function-terms energy-fn))
         (staple-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-staple)) terms)))
    staple-term))

(defun nonbond-term (energy-fn)
  (let* ((terms (chem:rigid-body-energy-function-terms energy-fn))
         (staple-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-nonbond)) terms)))
    staple-term))


(defun first-optimize (energy-function position &key trajectory)
  (chem:save-coordinates-from-vector energy-function position)
  (let* ((min (chem:make-minimizer energy-function))
         (nonbond-term (nonbond-term energy-function)))
    (chem:enable-print-intermediate-results min)
    (when trajectory
      (let* ((first-pseudo-aggregate (build-pseudo-aggregate energy-function position))
             (number-of-pseudo-atoms (chem:number-of-atoms first-pseudo-aggregate)))
        (chem:set-step-callback min
                                #'(lambda (pos &rest args)
                                    (let* ((pseudo-agg (build-pseudo-aggregate energy-function pos))
                                           (atom-positions (make-array (* 3 number-of-pseudo-atoms) :element-type 'single-float :adjustable nil)))
                                      (extract-coordinates-from-aggregate pseudo-agg atom-positions)
                                      (nglv::append-coordinates trajectory atom-positions))))))
    (cando:configure-minimizer min
                               :max-sd-steps 10000
                               :sd-tolerance 10000.0
                               :max-cg-steps 100
                               :cg-tolerance 1000.0
                               :max-tn-steps 0)
    (chem:disable nonbond-term)
    (restart-case
        (handler-bind
            ((ext:unix-signal-received (lambda (c) (print "Done") (invoke-restart 'skip-rest))))
          (chem:minimize min))
      (skip-rest ()
        :report "Skip rest of minimization"
        (chem:write-intermediate-results-to-energy-function min)
        (print "Skipping rest of minimization") ))
    (chem:load-coordinates-into-vector energy-function position)))


(eval-when (:compile-toplevel :load-toplevel :execute)
  (defparameter *interrupt-restart* nil))

(defmacro with-interruptable (() &rest body)
  (if *interrupt-restart*
      `(progn ,@body)
      `(let ((*interrupt-restart* t))
         (restart-case
             (handler-bind
                 ((ext:unix-signal-received (lambda (c) (print "Done") (invoke-restart 'skip-rest))))
               (progn
                 ,@body))
           (skip-rest ()
             :report "Skip rest of process"
             (chem:write-intermediate-results-to-energy-function min)
             (print "Skipping rest of process"))))))


(defun maybe-setup-minimizer-to-record-trajectory (minimizer energy-function position trajectory)
  (when trajectory
    (let* ((first-pseudo-aggregate (build-pseudo-aggregate energy-function position))
           (number-of-pseudo-atoms (chem:number-of-atoms first-pseudo-aggregate)))
      (chem:set-step-callback minimizer
                              #'(lambda (pos &rest args)
                                  (let* ((pseudo-agg (build-pseudo-aggregate energy-function pos))
                                         (atom-positions (make-array (* 3 number-of-pseudo-atoms) :element-type 'single-float :adjustable nil)))
                                    (extract-coordinates-from-aggregate pseudo-agg atom-positions)
                                    (nglv::append-coordinates trajectory atom-positions)))))))

(defun origami-optimize (energy-function position &key trajectory nonbond-on)
  (chem:save-coordinates-from-vector energy-function position)
  (let* ((min (chem:make-minimizer energy-function))
         (nonbond-term (nonbond-term energy-function)))
    (with-interruptable ()
      (chem:enable-print-intermediate-results min)
      (maybe-setup-minimizer-to-record-trajectory min energy-function position trajectory)
      (if nonbond-on
          (cando:configure-minimizer min
                                     :max-sd-steps 1000
                                     :sd-tolerance 10000000.0
                                     :max-cg-steps 1000
                                     :cg-tolerance 100.0
                                     :max-tn-steps 0)
          (cando:configure-minimizer min
                                     :max-sd-steps 1000
                                     :sd-tolerance 10000.0
                                     :max-cg-steps 1000
                                     :cg-tolerance 100.0
                                     :max-tn-steps 0))
      (if nonbond-on
          (chem:enable nonbond-term)
          (chem:disable nonbond-term))
      (cando:minimize-no-fail min))
    (chem:enable nonbond-term))
  (chem:load-coordinates-into-vector energy-function position))

(defun mc-criterion (delta delta-sum cycle temp)
  (exp (/ (/ (- delta) (/ delta-sum (float cycle)) temp))))

(defun monte-carlo (energyfn temp initial-pos cycles &key (step-size 1.0) verbose trajectory)
  (when trajectory
    (when (= (nglv::n-frames trajectory) 0)
      (let* ((pseudo-aggregate (origami::build-pseudo-aggregate energyfn initial-pos))
             (coordinates (make-array (* 3 (nglv::number-of-atoms trajectory)) :element-type 'single-float)))
        (extract-coordinates-from-aggregate pseudo-aggregate coordinates)
        (nglv::append-coordinates trajectory coordinates))))
  (let ((energy (chem:evaluate-energy energyfn initial-pos))
        (delta-sum 0.0)
        (keeps 0)
        (pos initial-pos)
        (rejects 0))
    (with-interruptable ()
      (loop for cycle from 1 to cycles
            for copy-pos = (copy-seq pos)
            for new-pos = (progn
                            (jostle-origami copy-pos step-size)
                            copy-pos)
            for new-energy = (chem:evaluate-energy energyfn new-pos)
            for kept = nil
            for delta = (- new-energy energy)
            do (incf delta-sum (abs delta))
            if (< delta 0.0)
              do (setf pos new-pos
                       energy new-energy
                       kept t)
                 (incf keeps)
                 (when verbose (format t "~a  ~a : ~a~%" cycle new-energy kept))
            else
              do (let ((random (random 1.0))
                       (mc-criterion (mc-criterion delta delta-sum cycle temp)))
                   (if (< mc-criterion random)
                       (progn
                         (setf pos new-pos
                               energy new-energy
                               kept t)
                         (incf keeps))
                       (incf rejects))
                   (when verbose
                     (format t "~a  ~a : ~a~%" cycle new-energy kept)
                     (format t "     mc(~a)   <  random(~a)~%" mc-criterion random)))
            when (and kept trajectory)
              do (let* ((pseudo-aggregate (origami::build-pseudo-aggregate energyfn pos))
                        (coordinates (make-array (length (nglv::get-coordinates trajectory 0)) :element-type 'single-float)))
                   (extract-coordinates-from-aggregate pseudo-aggregate coordinates)
                   (nglv::append-coordinates trajectory coordinates))))
    (when verbose
      (format t "monte carlo done ~a steps   energy: ~a    keeps: ~a     rejects: ~a~%"
              cycles energy keeps rejects))
    (map-into initial-pos #'identity pos)
    (values pos energy keeps rejects)))

  
