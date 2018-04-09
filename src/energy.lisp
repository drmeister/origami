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


(defun build-one-aggregate (strands energy-function)
  (let ((agg (chem:make-aggregate)))
    (loop for strand in strands
       for id = (order-id strand)
       for trans = (rigid-body-transform energy-function id)
       do (let ((m (build-one-molecule strand trans)))
            (chem:add-matter agg m)))
    agg))

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
      (loop for index from 0 below (chem:get-nvector-size energy) by 7
         do (progn
              (setf (elt start-pos (+ index 0)) (random-quat))
              (setf (elt start-pos (+ index 1)) (random-quat))
              (setf (elt start-pos (+ index 2)) (random-quat))
              (setf (elt start-pos (+ index 3)) (random-quat))
              (setf (elt start-pos (+ index 4)) (random-pos))
              (setf (elt start-pos (+ index 5)) (random-pos))
              (setf (elt start-pos (+ index 6)) (random-pos))))
      (chem:rigid-body-energy-function-normalize-position energy start-pos)
      (format t "start-pos -> ~a~%" start-pos)
      (chem:save-coordinates-from-vector energy start-pos))))

(defun move-rigid-body (energy-function id vec)
  (let ((pos (make-array (chem:get-nvector-size energy-function) :element-type 'double-float))
        (index (* id 7)))
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


(defun add-p-o-staple-term (ef id1 pos1 id2 pos2)
  (chem:energy-rigid-body-staple-add-term ef
                                          230.0 ; 230.0 ; OS-P from parm99.dat
                                          1.61 ; OS-p from parm99.dat
                                          (* 7 id1)
                                          pos1
                                          (* 7 id2)
                                          pos2))

(defun build-energy-rigid-body-staple (strands)
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
                                   (chem:get-position (chem:atom-with-name (residue p5-partner) :O3*)))))
       do (let ((p3-partner (forward-node p3)))
            (when (and p3-partner (<= my-order-id (order-id (strand p3-partner))))
              (add-p-o-staple-term ef
                                   my-order-id
                                   (chem:get-position (chem:atom-with-name (residue p3) :O3*))
                                   (order-id (strand p3-partner))
                                   (chem:get-position (chem:atom-with-name (residue p3-partner) :P)))))
       do (when p5hb
            (let ((p5hb-partner (forward-node p5hb)))
              (when (and p5hb-partner (<= my-order-id (order-id (strand p5hb-partner))))
                (add-p-o-staple-term ef
                                     my-order-id
                                     (chem:get-position (chem:atom-with-name (residue p5hb) :O3*))
                                     (order-id (strand p5hb-partner))
                                     (chem:get-position (chem:atom-with-name (residue p5hb-partner) :P))))))
       do (when p3hb
            (let ((p3hb-partner (backward-node p3hb)))
              (when (and p3hb-partner (<= my-order-id (order-id (strand p3hb-partner))))
                (add-p-o-staple-term ef
                                     my-order-id
                                     (chem:get-position (chem:atom-with-name (residue p3hb) :P))
                                     (order-id (strand p3hb-partner))
                                     (chem:get-position (chem:atom-with-name (residue p3hb-partner) :O3*)))))))
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


(defconstant +bdna-rise+ 3.32) ;; Angstroms

#| Mathematica code to plot enonbond function
Clear[rab];
a = { ax, ay, az};
b = {bx, by, bz};
rab[a_, b_] := Sqrt[(a - b).(a - b)]
enonbond[a_, b_, r0_, edep_] := 
 edep*((r0/rab[a, b])^12 - 2*(r0/rab[a, b])^6)
enonbond[{0, 0, 0}, {10, 0, 0}, 1.0, 0.1]
Plot[enonbond[{0, 0, 0}, {x, 0, 0}, 1.0, 0.01], {x, 0.8, 5}, 
PlotRange -> Full]
|#
(defun create-fast-nonbond-terms (nonbond-component strand index)
  (let* ((nodes (dna-strand-as-list-of-nodes strand))
         (num-nodes (length nodes))
         (total-height (* num-nodes +bdna-rise+))
         (start-height (* total-height -0.5)))
    (loop for node in (dna-strand-as-list-of-nodes strand)
       for z-height = (+ start-height (* index +bdna-rise+))
       for pos = (geom:vec 0.0 0.0 z-height)
       do (chem:energy-rigid-body-nonbond-set-term
           nonbond-component
           (prog1 index (incf index))
           (chem:make-atom :bogus :c)     ; no atom
           10.0      ; 10 Angstrom radius for a B-form DNA double helix
           0.01      ; Very little attraction
           0.0
           pos))
    index))


(defun build-fast-energy-rigid-body-nonbond (strands)
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
             for strand-nonbonds-count = (length (dna-strand-as-list-of-nodes strand))
             for running-count = strand-nonbonds-count then (+ running-count strand-nonbonds-count)
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
           do (setf start-index (create-fast-nonbond-terms nonbond-component strand start-index)))
        nonbond-component))))

(defun build-energy-function (origami)
  (let ((staple (build-energy-rigid-body-staple origami))
        (nonbonded (build-fast-energy-rigid-body-nonbond origami))
        (energy (build-rigid-body-energy-function origami)))
    (chem:rigid-body-energy-function-add-term energy staple)
    (chem:rigid-body-energy-function-add-term energy nonbonded)
    energy))
