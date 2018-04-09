(in-package :origami)

;;; ------------------------------------------------------------
;;;
;;;  double/single strand

(defclass dna-strand () 
  ((order-id :initarg :order-id :accessor order-id)
   (p5-end :initarg :p5-end :accessor p5-end)
   (p3-end :initarg :p3-end :accessor p3-end)))

(defclass double-stranded-dna (dna-strand) ())
(defclass single-stranded-dna (dna-strand) ())


(defun dna-strand-as-list-of-nodes (strand)
  (let ((p5 (p5-end strand))
        (p3 (p3-end strand)))
    (loop for cur = p5 then (forward-node cur)
       collect cur
       until (eq cur p3))))

(defun single-or-double-strand (vstrands)
  (let ((node-positions (make-hash-table :test #'eq)))
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for p-vec = (p-strand vstrand)
       for n-vec = (n-strand vstrand)
       do (locate-nodes node-positions p-vec)
       do (locate-nodes node-positions n-vec))   
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for p-vec = (p-strand vstrand)
       for n-vec = (n-strand vstrand)
       append (loop for index from 0 below (length p-vec)
                 for p-node = (elt p-vec index)
                 for n-node = (elt n-vec index)
                 when (or p-node n-node)
                 collect  (if p-node
                              (strand-chain p-node node-positions)
                              (let ((pair (strand-chain n-node node-positions)))
                                (list (second pair) (first pair))
                                pair))))))

(defun single-or-double-classified-strands (vstrands)
  (let ((pairs (single-or-double-strand vstrands)))
    (format t "pairs -> ~a~%" (mapcar (lambda (pair) (format nil "~a -> ~a~%" (identifier (first pair)) (identifier (second pair)))) pairs))
    (loop for pair in pairs
       for p5 = (first pair) 
       for p3 = (second pair)
       for order-id from 0
       collect (if (hbond-node p5)
                   (let ((strand (make-instance 'double-stranded-dna :p5-end p5 :p3-end p3 :order-id order-id)))
                     (setf (strand p5) strand
                           (strand (hbond-node p5)) strand
                           (strand p3) strand
                           (strand (hbond-node p3)) strand)
                     strand)
                   (let ((strand (make-instance 'single-stranded-dna :p5-end p5 :p3-end p3 :order-id order-id)))
                     (setf (strand p5) strand
                           (strand p3) strand)
                     strand)))))

(defun locate-nodes (hash-table vec)
  (loop for index from 0 below (length vec)
     for node = (elt vec index)
     when node
     do  (let ((vec-idx-pos (cons vec index)))
           (setf (gethash node hash-table) vec-idx-pos))))    

(defun strand-chain (node hash-table)
  (if (hbond-node node)
      (list (double-chain node #'backward-node #'forward-node hash-table) (double-chain node #'forward-node #'backward-node hash-table))
      (list (single-chain node #'backward-node hash-table) (single-chain node #'forward-node hash-table))))

(defun get-position (node hash-table)
  (let ((reference (gethash node hash-table)))
    (elt (car reference) (cdr reference))))

(defun (setf get-position) (new node hash-table)
  (let ((reference (gethash node hash-table)))
    (setf (elt (car reference) (cdr reference)) new)))

(defun single-chain (node accessor hash-table)
  (loop
     with i-node = node
     with first-node = node
     while (and (funcall accessor i-node)
                (not (eq (funcall accessor i-node) first-node))
                (not (hbond-node (funcall accessor i-node))))
     do (setf (get-position i-node hash-table) nil)
     do (setf i-node (funcall accessor i-node))
     finally (setf (get-position i-node hash-table) nil)
       (return i-node)))

(defun double-chain (node accessor1 accessor2 hash-table)
  (loop
     with first-node = node
     with i-node = node
     while (and (funcall accessor1 i-node)
                (not (eq (funcall accessor1 i-node) first-node))
                (hbond-node (funcall accessor1 i-node))
                (eq (hbond-node (funcall accessor1 i-node)) (funcall accessor2 (hbond-node i-node))))
     do (setf (get-position i-node hash-table) nil)
     do (setf (get-position (hbond-node i-node) hash-table) nil)
     do (setf i-node (funcall accessor1 i-node))        
     finally (setf (get-position i-node hash-table) nil)
       (setf (get-position (hbond-node i-node) hash-table) nil)
       (return i-node)))

;;; ------------------------------------------------------------
;;;
;;;  Geometry of cylinder


(defclass helix-constant ()
  ((radius :initform 1 :initarg :radius :accessor radius) ;Unit: nm
   (base-pair-angol :initform 190 :initarg :base-pair-angol :accessor base-pair-angol);Unit: degree
   (rotation-angol :initform 36 :initarg :rotation-angol :accessor rotation-angol);Unit: degree
   (rise :initform 3.2 :initarg :rise :accessor rise))) ;Unit: nm

(defparameter *bform-dna-transform* (make-instance 'helix-constant))
#|
(defun helix-coordinate-axis-points (constant start-index below-index step)
  (let ((initial-position (geom:vec 0.0 0.0 (* start-index (rise constant))))
        (points (make-array num)))
    (loop for i from start-index below below by step
         (
          (matrix (prodotto-matrici
                   (geom:make-m4-rotate-z (* 0.0174533 index (rotation-angol constant)))
                   (geom:make-m4-translate (list 0.0 0.0 (* index (rise constant)))))))
    (unless p-strand-sense
      (setf initial-position
            (geom:m*v (geom:make-m4-rotate-z
                       (* 0.0174533 (base-pair-angol constant)))
                       initial-position)))
    (geom:m*v matrix initial-position)))
|#
(defun helix-coordinate-base-model (constant index p-strand-sense)
  (let ((initial-position (geom:vec (radius constant) 0.0 0.0))
          (matrix (prodotto-matrici
                   (geom:make-m4-rotate-z (* 0.0174533 index (rotation-angol constant)))
                   (geom:make-m4-translate (list 0.0 0.0 (* index (rise constant)))))))
    (unless p-strand-sense
      (setf initial-position
            (geom:m*v (geom:make-m4-rotate-z
                       (* 0.0174533 (base-pair-angol constant)))
                       initial-position)))
    (geom:m*v matrix initial-position)))

(defun helix-transform (constant index)
  (geom:m*m (geom:make-m4-rotate-z (* 0.0174533 index (rotation-angol constant)))
            (geom:make-m4-translate (list 0.0 0.0 (* index (rise constant))))))


(defun extract-one-residue (aggregate residue-name)
  (let ((r0 (chem:content-at (chem:content-at aggregate 0) 0))
        (r1 (chem:content-at (chem:content-at aggregate 1) 0)))
    (if (eq residue-name (chem:get-name r0))
        r0
        r1)))

(defvar *bases* nil)

(defun load-bases (dir)
  "Load the DNA bases that are in the DNA origami"
  (let ((ht (make-hash-table))
        (cg (chem:load-pdb (namestring (make-pathname :name "C-G" :type "pdb" :defaults dir))))
        (gc (chem:load-pdb (namestring (make-pathname :name "G-C" :type "pdb" :defaults dir))))
        (at (chem:load-pdb (namestring (make-pathname :name "A-T" :type "pdb" :defaults dir))))
        (ta (chem:load-pdb (namestring (make-pathname :name "T-A" :type "pdb" :defaults dir)))))
    (let ((cg-c (extract-one-residue cg :C))
          (cg-g (extract-one-residue cg :G))
          (gc-g (extract-one-residue gc :G))
          (gc-c (extract-one-residue gc :C))
          (at-a (extract-one-residue at :A))
          (at-t (extract-one-residue at :T))
          (ta-t (extract-one-residue gc :T))
          (ta-a (extract-one-residue gc :A)))
      (setf (gethash :c ht) (cons cg-c cg-g))
      (setf (gethash :g ht) (cons gc-g gc-c))
      (setf (gethash :a ht) (cons at-a at-t))
      (setf (gethash :t ht) (cons ta-t ta-a)))
    (setf *bases* ht)))

(defun complementary-base-name (name)
  (cdr (assoc name '( ( :c . :g )
                     ( :g . :c )
                     ( :a . :t )
                     ( :t . :a )))))

(defun fill-residue (node transform &optional (bases *bases*))
  (let* ((base-name (base-name node))
         (residues (gethash base-name bases))
         (residue (chem:matter-copy (car residues))))
    (chem:apply-transform-to-atoms residue transform)
    (setf (residue node) residue)
    (when (hbond-node node)
      (let ((hbond-residue (chem:matter-copy (cdr residues))))
        (chem:apply-transform-to-atoms hbond-residue transform)
        (setf (residue (hbond-node node)) hbond-residue)))))

(defun fill-nodes-with-residues (strands &optional sequence)
  (loop for strand in strands
       do (loop for node in (dna-strand-as-list-of-nodes strand)
             for index from 0
             for transform = (helix-transform *bform-dna-transform* index)
             do (fill-residue node transform))))


(defun rigid-body-transform (energy-function id)
  (multiple-value-bind (a b c d x y z)
      (chem:rigid-body-energy-function-get-position energy-function id)
    (let ((m (geom:make-matrix nil)))
      (geom:set-from-quaternion m a b c d x y z)
      m)))


(defun load-origami (origami-pathname)
  (let* ((sin (open origami-pathname :direction :input))
         (json (json:decode-json sin))
         (origami (parse-json json)))
    origami))
