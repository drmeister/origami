
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
   (new-staple-vec :initarg :new-staple-vec :accessor new-staple-vec)
   (new-scaffold-vec :initarg :new-scaffold-vec :accessor new-scaffold-vec)
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
       do (intra-helix-connect (staple-vec vstrand) (scaffold-vec vstrand));to be postpone after skip-loop 
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
	    (setf (new-staple-vec vstrand) a
		  (new-scaffold-vec vstrand) b)))   
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
    (loop for step in step-direction
	with pos = 0 and neg = 0
	when step
	do (cond ((plusp step) (incf pos))
		 ((minusp step) (incf neg))
		 (T (error "Do we find a basis connected with itself? ~a" step-direction)))
       finally (return
		 (cond ((> pos neg) 1)
		       ((< pos neg) -1)
		       (T (error "arrow direction in the vector is neither forward nor backward ~a" step-direction)))))))

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
	     (new-scaf-vec (make-array new-vec-length))
	     (new-stap-vec (make-array new-vec-length)))
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
		       (cond
			 ((and (= stap-direction 1)
			       (= scaf-direction -1))
					;stap->f
					;scaf->b
;			  (format t "stap->f scaf->b~%")
			  (and (elt staple-vec s-l-cursor) (loop-procedure num :ST staple-vec s-l-cursor new-stap-vec dest-cursor loop-vec #'(setf backward-node) #'backward-node #'(setf forward-node)))
			  (and (elt scaffold-vec s-l-cursor) (loop-procedure num :SC scaffold-vec s-l-cursor new-scaf-vec dest-cursor loop-vec #'(setf forward-node) #'forward-node #'(setf backward-node))))
			 ((and (= stap-direction -1)
			       (= scaf-direction 1))
					;stap->b
					;scaf->f
;			  (format t "stap->b scaf->f~%")
			  (and (elt staple-vec s-l-cursor) (loop-procedure num :ST staple-vec s-l-cursor new-stap-vec dest-cursor loop-vec #'(setf forward-node) #'forward-node #'(setf backward-node)))
			  (and (elt scaffold-vec s-l-cursor) (loop-procedure num :SC scaffold-vec s-l-cursor new-scaf-vec dest-cursor loop-vec #'(setf backward-node) #'backward-node #'(setf forward-node))))
			 (t (error "What do I do with stap-direction = ~a and scaf-direction = ~a~%"
				   stap-direction scaf-direction)))
		       (let ((staple-node (elt new-stap-vec dest-cursor))
			     (scaffold-node (elt new-scaf-vec dest-cursor)))
			 (when (and staple-node scaffold-node)
			   (setf (hbond-node staple-node) scaffold-node)
			   (setf (hbond-node scaffold-node) staple-node)))
		       (decf (elt loop-vec s-l-cursor))
		       (incf dest-cursor))
		      ((and (= (elt skip-vec s-l-cursor) 0)
			    (= (elt loop-vec s-l-cursor) 0))
		       ;;do a copy
;		       (format t "copy~%")
		       (setf (elt new-stap-vec dest-cursor) (elt staple-vec s-l-cursor))
		       (setf (elt new-scaf-vec dest-cursor) (elt scaffold-vec s-l-cursor))
		       (incf dest-cursor)
		       (incf s-l-cursor))
		      (t (error "An impossible step was encountered")))
		    (incf s-l-cursor)
		    )))
	(values new-stap-vec new-scaf-vec)
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

(defun draw-result-graph (dest strands)
  (let ((dest (if (streamp dest)
		  dest
		  (open dest :direction :output :if-exists :supersede))))
;    (format t "1")
    (format dest "digraph G {~%")
    (loop for strand being the hash-values in strands using (hash-key num)
       for new-staple-vec = (new-staple-vec strand)
       for new-scaffold-vec = (new-scaffold-vec strand)
;       do (format t "2")
       do (draw-double-strand dest num new-scaffold-vec new-staple-vec))
    (loop for strand being the hash-values in strands using (hash-key num)
       for new-staple-vec = (new-staple-vec strand)
       for new-scaffold-vec = (new-scaffold-vec strand)
       do (loop for node across new-scaffold-vec
	     when node
	     do (draw-connection dest node))
       do (loop for node across new-staple-vec
	     when node
	     do (draw-connection dest node))
	 )
    (format dest "}~%")))

#| testing code

	(defparameter result (parse-json *j*))

	result

	(/= -1 -1 -2 -3 0)

	(scaffold-vec (gethash 0 *origami*))

        (skip-loop (skip-json (gethash 0 *origami*)) (loop-json (gethash 0 *origami*)) (staple-vec (gethash 0 *origami*)) (scaffold-vec (gethash 0 *origami*)))


	|#

