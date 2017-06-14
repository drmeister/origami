
(ql:quickload "cl-json")


(in-package :cl-jupyter-user)

(defun LIST-OF-VSTRANDS-FROM-JSON (json)
  (cdr (assoc :vstrands json)))


(defclass vstrand ()
  ((staple-json :initarg :staple-json :reader staple-json)
   (scaffold-json :initarg :scaffold-json :reader scaffold-json)
   (skip-json :initarg :skip-json :reader skip-json)
   (loop-json :initarg :loop-json :reader loop-json)
   (scaf-loop-json :initarg :scaf-loop-json :reader scaf-loop-json)
   (stap-loop-json :initarg :stap-loop-json :reader stap-loop-json)
   (staple-vec :initarg :staple-vec :reader staple-vec)
   (scaffold-vec :initarg :scaffold-vec :reader scaffold-vec)
   (p-strand :initarg :p-strand :accessor p-strand)
   (n-strand :initarg :n-strand :accessor n-strand)
   (row :initarg :row :reader row)
   (col :initarg :column :reader column)
   ))

  ;check if not ("format":"3.0")
(defun parse-json (json)
  (let ((vstrands (make-hash-table :test #'eql)))
    (loop for json-vstrand-info in (LIST-OF-VSTRANDS-FROM-JSON json)
       for num = (cdr (assoc :num json-vstrand-info))
       for col = (cdr (assoc :col json-vstrand-info))
       for row = (cdr (assoc :row json-vstrand-info))
       for staple-json = (cdr (assoc :stap json-vstrand-info))
       for scaffold-json = (cdr (assoc :scaf json-vstrand-info))
       for skip = (cdr (assoc :skip json-vstrand-info))
       for loop = (cdr (assoc :loop json-vstrand-info))
       for scaf-loop = (cdr (assoc :scafLoop json-vstrand-info))
       for stap-loop = (cdr (assoc :stapLoop json-vstrand-info))
       for vstrand = (progn
		       (or (= (length staple-json) (length scaffold-json))
			   (error "In vstrand ~a the staple length ~a and the scaffold length ~a are not the same - and they must be"
				  num
				  (length staple-json)
				  (length scaffold-json)))
		       (make-instance 'vstrand
				      :staple-json staple-json
				      :scaffold-json scaffold-json
				      :skip-json skip
				      :loop-json loop
				      :scaf-loop-json scaf-loop
				      :stap-loop-json stap-loop
				      :staple-vec (make-array (length staple-json))
				      :scaffold-vec (make-array (length scaffold-json))
				      :row row
				      :column col))
       do (format t "vstrand -> ~a~%" vstrand)
       do (setf (gethash num vstrands) vstrand)
       do (BUILD-NODE staple-json (staple-vec vstrand) num :sT)
       do (BUILD-NODE scaffold-json (scaffold-vec vstrand) num :sC)
       do (intra-helix-connect (staple-vec vstrand) (scaffold-vec vstrand))
	 )
    ;; connect the nodes
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for staple-json = (staple-json vstrand)
       for staple-vec = (staple-vec vstrand)
       for scaffold-json = (scaffold-json vstrand)
       for scaffold-vec = (scaffold-vec vstrand)
       do (CONNECT-EVERYTHING vstrands num staple-json staple-vec #'staple-vec)
       do (CONNECT-EVERYTHING vstrands num scaffold-json scaffold-vec #'scaffold-vec))
    ;; Deal with loop
    ;; (skip-loop skip loop staple-vec scaffold-vec) -> return new staple and scaffold
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for loop-json = (loop-json vstrand)
       for staple-vec = (staple-vec vstrand)
       for skip-json = (skip-json vstrand)
       for scaffold-vec = (scaffold-vec vstrand)
       do (format t "duplex num ~a~%" num)
       do (multiple-value-bind (a b)
	      (skip-loop num skip-json loop-json staple-vec scaffold-vec)
	    (setf (p-strand vstrand) a
		  (n-strand vstrand) b)))   
    vstrands))

    
(defun BUILD-NODE (strand-json vec num strand-name)
  (loop for index from 0 below (length vec)
     for node-json in strand-json
     unless (apply #'= -1 node-json)
     do (setf (elt vec index) (make-instance 'node :name (list :j num strand-name index)))))

(defun INTRA-HELIX-CONNECT (staple-vec scaffold-vec)
  (loop for index from 0 below (length staple-vec)
     for staple-node = (elt staple-vec index)
     for scaffold-node = (elt scaffold-vec index)
     do (when (and staple-node scaffold-node)
	  (setf (hbond-node staple-node) scaffold-node)
	  (setf (hbond-node scaffold-node) staple-node))))
    
(defun connect-everything (vstrands num strand-json strand-vec accessor)
  (loop for json-info in strand-json
     for index from 0 below (length strand-json)
     for node = (elt strand-vec index)
     when node
     do (let* ((fwd-vstrand (first json-info))
	       (fwd-pos     (second json-info))
	       (bkd-vstrand (third json-info))
	       (bkd-pos     (fourth json-info)))
	  (when (/= fwd-vstrand -1)
	    (let ((forward-node (lookup-node vstrands fwd-vstrand fwd-pos accessor)))
	      (assert forward-node)
	      (setf (forward-node node) forward-node)))
	  (when (/= bkd-vstrand -1)
	    (let ((backward-node (lookup-node vstrands bkd-vstrand bkd-pos accessor)))
	      (assert backward-node)
	      (setf (backward-node node) backward-node))))))
	      
(defclass node ()
  ((hbond-node :initform nil :initarg :hbond-node :accessor hbond-node)
   (t-q-bond-node :initform nil :initarg :t-q-bond-node :accessor t-q-bond-node)
   (forward-node :initform nil :initarg :forward-node :accessor forward-node)
   (backward-node :initform nil :initarg :backward-node :accessor backward-node)
   (name :initarg :name :reader name)))

(defun lookup-node (vstrands vstrand-num pos accessor)
  (let* ((vstrand (gethash vstrand-num vstrands))
	 (strand (funcall accessor vstrand)))
    (elt strand pos)))


(defun arrow-direction (vec)
  (let ((step-direction (loop for index from 0 below (length vec)
			   for node = (elt vec index)
			   when node
			   collect (let ((fwd (forward-node node)))
				     (when (and fwd (position fwd vec))
				       (- (position fwd vec) index))))))
#|    (cond
      ((every (lambda (x) (if x (plusp x) T)) step-direction) 1)
      ((every (lambda (x) (if x (minusp x) T)) step-direction) -1)
      (t (error "arrow direction in the vector is neither forward nor backward ~a" step-direction)))))|#
					;(if (every nil step-direction) (error "empty vector")
    (if step-direction 
	(loop for step in step-direction
	   with pos = 0 and neg = 0
	   when step
	   do (cond ((plusp step) (incf pos))
		    ((minusp step) (incf neg))
		    (T (error "Do we find a basis connected with itself? ~a" step-direction)))
	   finally (return (cond ((> pos neg) 1)
				 ((< pos neg) -1)
				 (T (error "arrow direction in the vector is neither forward nor backward ~a" step-direction)))))
	nil)))

(defun skip-loop (num skip-json loop-json staple-vec scaffold-vec)
  (flet ((skip-procedure (x)
	   (let ((xf (forward-node x))
		 (xb (backward-node x)))
	     (when xb (setf (forward-node xb) xf))
	     (when xf (setf (backward-node xf) xb))
	     (setf x NIL)))
	 (loop-procedure (num strand-name old-vec s-l-cursor new-vec dest-cursor loop-vec setf-accessor-l accessor-l setf-accessor-r)
	    (setf (elt new-vec dest-cursor) (make-instance 'node :name (list :A num strand-name s-l-cursor :I (elt loop-vec s-l-cursor))))
	    (let ((New-node (elt new-vec dest-cursor))
		  (Doubled-node (elt old-vec s-l-cursor))
		  (l-d-node (funcall accessor-l (elt old-vec s-l-cursor))))
	      (when Doubled-node (funcall setf-accessor-l New-node Doubled-node))
	      (funcall setf-accessor-r  Doubled-node New-node)
	      (funcall setf-accessor-l l-d-node New-node)
	      (when l-d-node (funcall setf-accessor-r New-node l-d-node)))))
    ;;must check if 4 input have the same lenght
    (let ((skip-vec (coerce skip-json 'vector))
	  (loop-vec (coerce loop-json 'vector))
	  (stap-direction (arrow-direction staple-vec))
	  (scaf-direction (arrow-direction scaffold-vec)))
      (format t "Starting function~%")
      (let* ((old-vec-min-length
	      (loop for index from 0 below (length scaffold-vec)
		 when (or (elt scaffold-vec index) (elt staple-vec index))
		 count 1))
	     (new-vec-length (+ old-vec-min-length (reduce #'+ skip-json) (reduce #'+ loop-json)))
	     (p-strand (make-array new-vec-length));p-strand
	     (n-strand (make-array new-vec-length))
	     (condition (cond
			  ((and (not stap-direction)
				(not scaf-direction))
			   t)
			  ((and (or (not stap-direction) (= stap-direction 1))
				(or (not scaf-direction) (= scaf-direction -1)))
			   t)
			  ((and (or (not stap-direction) (= stap-direction -1))
				(or (not scaf-direction) (= scaf-direction 1)))
			   nil)
			  (t (error "What do I do with stap-direction = ~a and scaf-direction = ~a~%" stap-direction scaf-direction))))
	     (old-p-strand (if condition staple-vec scaffold-vec))
	     (old-n-strand (if condition scaffold-vec staple-vec)))
	(format t "Starting loop~%")
	(loop with s-l-cursor = 0
	   with dest-cursor = 0
	   with length-skip = (length skip-vec)
	   until (= s-l-cursor length-skip)
	   do (progn
;		(format t "s-l-cursor -> ~a~%" s-l-cursor)
		(if (or (elt scaffold-vec s-l-cursor) (elt staple-vec s-l-cursor))
		    (cond
		      ((and (= (elt skip-vec s-l-cursor) -1)
			    (= (elt loop-vec s-l-cursor) 0))
		       ;;do a skip procedure
		       (format t "hit a skip~%")
		       (skip-procedure (elt staple-vec s-l-cursor))
		       (skip-procedure (elt scaffold-vec s-l-cursor))
		       (incf s-l-cursor)
		       )
		      ((and (= (elt skip-vec s-l-cursor) 0)
			    (> (elt loop-vec s-l-cursor) 0))
		       ;;do an insertion
		       (format t "Insertion number -> ~a~%" (elt loop-vec s-l-cursor))
		       (and (elt staple-vec s-l-cursor)
			    (loop-procedure num :ST old-p-strand s-l-cursor p-strand dest-cursor loop-vec #'(setf backward-node) #'backward-node #'(setf forward-node))); p-strand
		       (and (elt scaffold-vec s-l-cursor)
			    (loop-procedure num :SC old-n-strand s-l-cursor n-strand dest-cursor loop-vec #'(setf forward-node) #'forward-node #'(setf backward-node))); n-strand
		      (let ((p-node (elt p-strand dest-cursor)); p-node p-strand
			     (n-node (elt n-strand dest-cursor))); n-node n-strand
			 (when (and p-node n-node); p-node n-node
			   (setf (hbond-node p-node) n-node); p-node n-node
			   (setf (hbond-node n-node) p-node))); n-node p-node 
		       (decf (elt loop-vec s-l-cursor))
		       (incf dest-cursor))
		      ((and (= (elt skip-vec s-l-cursor) 0)
			    (= (elt loop-vec s-l-cursor) 0))
		       ;;do a copy
;		       (format t "copy~%")
		       (setf (elt p-strand dest-cursor) (elt old-p-strand s-l-cursor)); n or p
		       (setf (elt n-strand dest-cursor) (elt old-n-strand s-l-cursor));n or p
		       (incf dest-cursor)
		       (incf s-l-cursor))
		      (t (error "An impossible step was encountered")))
		    (incf s-l-cursor)
		    )))
	(values p-strand n-strand); n-strand p-strand
;	(list new-stap-vec new-scaf-vec)
;	new-stap-vec
	))))




;;; ------------------------------------------------------------
;;;
;;;  Graphviz generated

(defun safe-draw-link (dest node next-node color &optional link)
  ;(let ((*print-pretty* nil)))
  (when (and node next-node)
    (format dest "\"~a\" -> \"~a\" [color=~a];~%" (name node) (name next-node) color)))

(defun draw-connection (dest node)
  (safe-draw-link dest node (forward-node node) "violet")
  (safe-draw-link dest node (backward-node node) "blue")
  (safe-draw-link dest node (hbond-node node) "orange"))

(defun draw-node (dest node)
;  (let ((*print-pretty* nil)))
  (format dest "\"(~{~a~^ ~})\"->" (name node)))

(defun draw-strand (dest num name vec)
  (format dest "subgraph cluster~a_~a {~%" (string name) num )
  (format dest "     edge [style=\"invis\"];~%")
  (format dest "     label = \"~a\";~%" (string name))
  (format dest "     color = blue;~%")
  (format dest "     \"start~a_~a\"  [style=\"invis\"];~%" (string name) num)
  (format dest "     \"stop~a_~a\"  [style=\"invis\"];~%" (string name) num)
  (format dest "     \"start~a_~a\"->" (string name) num)
  ;(direction (arrow-direction vec))
  (loop for node across vec
     ;for direction =(arrow-direction vec); (if (= 1 (arrow-direction vec)) (#'forward-node) (#'backward-node))
     when node
     do (draw-node dest node))
  (format dest "\"stop~a_~a\"}~%" (string name) num))
  
(defun draw-double-strand (dest num scaffold staple)
  (format dest "subgraph cluster_~a {~%" num)
  (format dest "label = \"double-strand-~a\";~%" num)
  (draw-strand dest num :scaffold scaffold)
  (draw-strand dest num :staple staple)
  (format dest "}~%"))


(defun draw-graph (dest strands)
  (let ((dest (if (streamp dest)
		  dest
		  (open dest :direction :output :if-exists :supersede))))
    (format dest "digraph G {~%")
    (loop for strand being the hash-values in strands using (hash-key num)
       for staple-vec = (staple-vec strand)
       for scaffold-vec = (scaffold-vec strand)
       do (draw-double-strand dest num scaffold-vec staple-vec))
    (loop for strand being the hash-values in strands using (hash-key num)
       for staple-vec = (staple-vec strand)
       for scaffold-vec = (scaffold-vec strand)
       do (loop for node across scaffold-vec
	     when node
	     do (draw-connection dest node))
       do (loop for node across staple-vec
	     when node
	     do (draw-connection dest node))
	 )
    (format dest "}~%")))

(defun draw-result-double-strand (dest num p-strand n-strand)
  (format dest "subgraph cluster_~a {~%" num)
  (format dest "label = \"double-strand-~a\";~%" num)
  (draw-strand dest num :p p-strand)
  (draw-strand dest num :n n-strand)
  (format dest "}~%"))

(defun draw-result-graph (dest strands)
  (let ((dest (if (streamp dest)
		  dest
		  (open dest :direction :output :if-exists :supersede))))
;    (format t "1")
    (format dest "digraph G {~%")
    (loop for strand being the hash-values in strands using (hash-key num)
       for p-strand = (p-strand strand)
       for n-strand = (n-strand strand)
;       do (format t "2")
       do (draw-double-strand dest num p-strand n-strand))
    (loop for strand being the hash-values in strands using (hash-key num)
       for p-strand = (p-strand strand)
       for n-strand = (n-strand strand)
       do (loop for node across p-strand
	     when node
	     do (draw-connection dest node))
       do (loop for node across n-strand
	     when node
	     do (draw-connection dest node))
	 )
    (format dest "}~%")))

;;; ------------------------------------------------------------
;;;
;;;  double/single strand

(defun single-or-double-strand (vstrands)
  (let ((node-positions (make-hash-table :test #'eq)))
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for st-vec = (p-strand vstrand)
       for sc-vec = (n-strand vstrand)
       do (locate-nodes node-positions st-vec)
       do (locate-nodes node-positions sc-vec))	 
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
     for st-vec = (p-strand vstrand)
     for sc-vec = (n-strand vstrand)
     append (loop for index from 0 below (length st-vec)
	   for st-node = (elt st-vec index)
	   for sc-node = (elt sc-vec index)
	   when (or st-node sc-node)
	   collect  (if st-node
			(strand-chain st-node node-positions)
			(strand-chain sc-node node-positions))))))

(defun locate-nodes (hash-table vec)
  (loop for index from 0 below (length vec)
     for node = (elt vec index)
     when node
     do  (let ((vec-idx-pos (cons vec index)))
	   (setf (gethash node hash-table) vec-idx-pos))))    

(defun strand-chain (node hash-table)
  (if (hbond-node node)
      (list (double-chain node #'forward-node #'backward-node hash-table) (double-chain node #'backward-node #'forward-node hash-table))
      (list (single-chain node #'forward-node hash-table) (single-chain node #'backward-node hash-table))))

(defun get-position (node hash-table)
  (let ((reference (gethash node hash-table)))
    (elt (car reference) (cdr reference))))

(defun (setf get-position) (new node hash-table)
  (let ((reference (gethash node hash-table)))
    (setf (elt (car reference) (cdr reference)) new)))

(defun single-chain (node accessor hash-table)
  (loop with i-node = node
     while (and (funcall accessor i-node) (not (hbond-node (funcall accessor i-node))))
     do (setf (get-position i-node hash-table) nil)
     do (setf i-node (funcall accessor i-node))
     finally (setf (get-position i-node hash-table) nil)
       (return i-node)))


(defun double-chain (node accessor1 accessor2 hash-table)
  (loop with i-node = node
     while (and (funcall accessor1 i-node) (hbond-node (funcall accessor1 i-node))
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
   (rotation-angol :initform 34.3 :initarg :rotation-angol :accessor rotation-angol);Unit: degree
   (rise :initform 3.32 :initarg :rise :accessor rise))) ;Unit: nm

(defun helix-coordinate (constant index p-strand-sense)
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



#| testing code

	(defparameter result (parse-json *j*))

	result

	(/= -1 -1 -2 -3 0)

	(scaffold-vec (gethash 0 *origami*))

        (skip-loop (skip-json (gethash 0 *origami*)) (loop-json (gethash 0 *origami*)) (staple-vec (gethash 0 *origami*)) (scaffold-vec (gethash 0 *origami*)))


	|#

