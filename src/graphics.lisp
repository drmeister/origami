(in-package :origami)


;;; Generate an nglview shape for the nonbond term

#+(or)(defun one-cylinder (centers num)
  (let ((p0 (first centers))
        (pn (car (last centers))))
    (vector "cylinder" 
            (vector (float (first p0) 1.0s0)
                    (float (second p0) 1.0s0)
                    (float (third p0) 1.0s0))
            (vector (float (first pn) 1.0s0)
                    (float (second pn) 1.0s0)
                    (float (third pn) 1.0s0))
            (vector 1.0 0.0 0.0)
            3.0)))

(defun one-cylinder (positions)
  (let* ((sphere-list (loop for pos in positions
                            for num from 0
                            collect (vector "sphere" 
                                            (vector (float (first pos) 1.0s0)
                                                    (float (second pos) 1.0s0)
                                                    (float (third pos) 1.0s0))
                                            (if (= num 0)
                                                (vector 0.0 1.0 0.0)
                                                (vector 1.0 0.0 0.0))
                                            10.0)))
         (sphere-vec (coerce sphere-list 'vector)))
    sphere-vec))


(defun nonbond-term-points (energy-fn position)
  (let* ((terms (chem:rigid-body-energy-function-terms energy-fn))
         (nonbond-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-nonbond))
                                terms))
         (parts (chem:parts-as-list nonbond-term position))
         (cylinders (mapcar (lambda (one)
                              (mapcar (lambda (x)
                                        (list (float (first x) 1.0s0)
                                              (float (second x) 1.0s0)
                                              (float (third x) 1.0s0)))
                                      one))
                            parts)))
    cylinders))

(defun staple-term-parts (energy-fn position)
  (let* ((terms (chem:rigid-body-energy-function-terms energy-fn))
         (staple-term (find-if (lambda (x) (typep x 'chem:energy-rigid-body-staple)) terms)))
    (when staple-term
      (chem:parts-as-list staple-term position))))

(defun add-spheres (widget positions)
  (let ((spheres (coerce (loop for pos in positions
                               for num from 0
                               collect (vector "sphere" 
                                               (vector (float (first pos) 1.0s0)
                                                       (float (second pos) 1.0s0)
                                                       (float (third pos) 1.0s0))
                                               (if (= num 0)
                                                   (vector 0.0 1.0 0.0)
                                                   (vector 1.0 0.0 0.0))
                                               5.0))
                         'vector)))
         (cando::safe-add-shape widget spheres :name (format nil "cylinder~a" num))
         spheres))

(defun add-shape-for-nonbond-term-one-helix (widget energy-fn position num)
  (let* ((points (nonbond-term-points energy-fn position))
         (cylinder (elt points num)))
    (cando::safe-add-shape widget (one-cylinder cylinder) :name (format nil "cyl~a" num))))

(defun add-shape-for-nonbond-term (widget energy-fn position)
  (let ((cylinders (nonbond-term-points energy-fn position)))
    (loop for num from 0 below (length cylinders)
          for cyl in cylinders
          for cyl-shape = (one-cylinder cyl)
          do (cando::safe-add-shape widget cyl-shape :name (format nil "cyl~a" num)))))

(defun create-shape-for-staples (energy-function position)
  (let ((parts (staple-term-parts energy-function position)))
    (coerce (loop for part in parts
                  for x1 = (float (elt part 0) 1.0s0)
                  for y1 = (float (elt part 1) 1.0s0)
                  for z1 = (float (elt part 2) 1.0s0)
                  for x2 = (float (elt part 3) 1.0s0)
                  for y2 = (float (elt part 4) 1.0s0)
                  for z2 = (float (elt part 5) 1.0s0)
                  for r0 = (float (elt part 6) 1.0s0)
                  collect (vector "cylinder" (vector x1 y1 z1)
                                  (vector x2 y2 z2)
                                  (vector 0.0 1.0 1.0)
                                  0.5))
            'vector)))


(defun add-shape-for-staples (widget energy pos)
  (let ((shape (create-shape-for-staples energy pos)))
    (cando::safe-add-shape widget shape :name "staples")))
  
(defun add-shape-for-centers (widget pos)
  (let* ((spheres (loop for num from 0 below (/ (length pos) 7)
                       for x = (elt pos (+ (* num 7) 4))
                       for y = (elt pos (+ (* num 7) 5))
                       for z = (elt pos (+ (* num 7) 6))
                       collect (vector "sphere"
                                       (vector (float x 1.0s0)
                                               (float y 1.0s0)
                                               (float z 1.0s0))
                                       (vector 0.0 0.0 1.0)
                                       11.0)))
         (sphere-vec (coerce spheres 'vector)))
    (cando::safe-add-shape widget sphere-vec :name "centroid")))
                       

(defun add-shapes-for-energy (widget energy position)
  (add-shape-for-nonbond-term widget energy position)
  (add-shape-for-centers widget position)
  (add-shape-for-staples widget energy position))
