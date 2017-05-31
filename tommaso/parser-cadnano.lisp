
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
					;       do (BUILD-NODE (staple-vec vstrand) staple-json)
       do (BUILD-NODE staple-json (staple-vec vstrand) num :staple)
					;       do (BUILD-NODE (scaffold-vec vstrand) scaffold-json)
       do (BUILD-NODE scaffold-json (scaffold-vec vstrand) num :scaffold)
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
    vstrands))
    
(defun BUILD-NODE (strand-json vec num strand-name)
  (loop for index from 0 below (length vec)
     for node-json in strand-json
     unless (apply #'= -1 node-json)
     do (setf (elt vec index) (make-instance 'node :name (list :json num strand-name index)))))

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
				     (when fwd
				       (- (position fwd vec) index))))))
    (cond
      ((every (lambda (x) (if x (plusp x) T)) step-direction) 1)
      ((every (lambda (x) (if x (minusp x) T)) step-direction) -1)
      (t (error "arrow direction in the vector is neither forward nor backward")))))



(defun skip-loop (skip-json loop-json staple-vec scaffold-vec)
  (flet ((skip-procedure (x)
	   (let ((xf (forward x))
		 (xb (backward x)))
	     (setf (forward xb) xf)
	     (setf (backward xf) xb)
	     (setf x NIL)))
	 (loop-procedure (old-vec s-l-cursor new-vec dest-cursor loop-vec setf-accessor-l setf-accessor-r)
	    (setf (elt new-vec dest-cursor) (make-instance 'node :name (cons :loop (elt loop-vec s-l-cursor))))
	    (let ((New-node (elt new-vec dest-cursor))
		  (Doubled-node (elt old-vec s-l-cursor))
		  (l-d-node (accessor-l (elt old-vec s-l-cursor))))
	      (funcall setf-accessor-l New-node Doubled-node)
	      (funcall setf-accessor-r  Doubled-node New-node)
	      (funcall setf-accessor-l l-d-node New-node)
	      (funcall setf-accessor-r New-node l-d-node))))
    ;;must check if 4 input have the same lenght
    (let ((skip-vec (coerce skip-json 'vector))
	  (loop-vec (coerce loop-json 'vector))
	  (stap-direction (arrow-direction staple-vec))
	  (scaf-direction (arrow-direction scaffold-vec)))
      (let* ((old-vec-min-length
	      (loop for index from 0 below (length scaffold-vec)
		 when (or (elt scaffold-vec index) (elt staple-vec index))
		 count))
	     (new-vec-length (apply #'+ old-vec-min-length (apply #'+ skip-vec) loop-vec))
	     (new-scaf-vec (make-array new-vec-length))
	     (new-stap-vec (make-array new-vec-length)))
	#|
	(cond ((and (= stap-direction 1)
	(= scaf-direction -1))
					;stap->f ; ; ; ;
					;scaf->b ; ; ; ;
	)
	((and (= stap-direction -1)
	(= scaf-direction 1))
					;stap->b ; ; ; ;
					;scaf->f ; ; ; ;
	)
	(t (error "What do I do with stap-direction = ~a and scaf-direction = ~a~%"
	stap-direction scaf-direction)))
	|#
    (loop with s-l-cursor = 0
       with dest-cursor = 0
       with length-skip = (length skip-vec)
       until (= s-l-cursor length-skip)
       do (cond
	    ((and (= (elt skip-vec s-l-cursor) -1)
		  (= (elt loop-vec s-l-cursor) 0))
	     ;;do a skip procedure
	     (skip-procedure (elt staple-vec s-l-cursor))
	     (skip-procedure (elt scaffold-vec s-l-cursor))
	     (incf s-l-cursor)
	     )
	    ((and (> (elt skip-vec s-l-cursor) 0)
		  (= (elt loop-vec s-l-cursor) 0))
	     ;;do an insertion
	     (cond
	       ((and (= stap-direction 1)
			 (= scaf-direction -1))
					;stap->f
					;scaf->b
		    (loop-procedure stap-vec s-l-cursor new-stap-vec dest-cursor #'(setf backward-node) #'(setf forward-node))
		    (loop-procedure scaffold-vec s-l-cursor new-stap-vec dest-cursor #'(setf forward-node) #'(setf backward-node)))
	       ((and (= stap-direction -1)
			 (= scaf-direction 1))
					;stap->b
					;scaf->f
		    (loop-procedure stap-vec s-l-cursor new-stap-vec dest-cursor #'(setf forward-node) #'(setf backward-node))
		    (loop-procedure scaffold-vec s-l-cursor new-stap-vec dest-cursor #'(setf backward-node) #'(setf forward-node)))
	       (t (error "What do I do with stap-direction = ~a and scaf-direction = ~a~%"
			     stap-direction scaf-direction)))
	     #|
	(let ((New-node (elt new-f-vec dest-cursor))
	     (Doubled-node (elt f-vec s-l-cursor))
	     (b-d-node (backward (elt f-vec s-l-cursor))))
	     (setf New-node (make-instance 'node :name (cons :loop (elt loop s-l-cursor))))
	     (setf (backward Doubled-node) New-node)
	     (setf (forward New-node) Doubled-node)
	     (setf (backward New-node) b-d-node)
	     (setf (forward b-d-node) New-node)
	     )
	     ;;backward arrow		; ; ;
	(let ((New-node (elt new-b-vec dest-cursor))
	     (Doubled-node (elt b-vec s-l-cursor))
	     (f-d-node (forward (elt b-vec s-l-cursor))))
	     (setf New-node (make-instance 'node :name (cons :loop (elt loop s-l-cursor))))
	     (setf (forward Doubled-node) New-node)
	     (setf (backward New-node) Doubled-node)		 
	     (setf (forward New-node) f-d-node)
	     (setf (backward b-d-node) New-node)
	     )
	|#
	     (decf (elt loop-vec s-l-cursor))
	     (incf dest-cursor))
	    ((and (= (elt skip-vec s-l-cursor) 0)
		  (= (elt loop-vec s-l-cursor) 0))
	     ;;do a copy
	     (setf (elt new-stap-vec dest-cursor) (elt staple-vec s-l-cursor))
	     (setf (elt new-scaf-vec dest-cursor) (elt scaffold-vec s-l-cursor))
	     (incf dest-cursor)
	     (incf s-l-cursor))
	    (t (error "An impossible step was encountered"))))))))




;;; ------------------------------------------------------------
;;;
;;;  Graphviz generated
(defun safe-draw-link (dest node next-node &optional link)
  (let ((*print-pretty* nil))
  (when (and node next-node)
    (format dest "\"~a\" -> \"~a\";~%" (name node) (name next-node)))))

(defun draw-node (dest node)
  (safe-draw-link dest node (forward-node node))
  (safe-draw-link dest node (backward-node node))
  (safe-draw-link dest node (hbond-node node)))
 
(defun draw-strand (dest num name vec)
  (format dest "subgraph strandv~a_~a {~%" (string name) num )
  (format dest "     label = \"~a\";~%" (string name))
  (format dest "     color = blue;~%")
  (loop for node across vec
     when node
     do (draw-node dest node))
  (format dest "}~%"))
  
(defun draw-double-strand (dest num scaffold staple)
  (format dest "subgraph doubleVstrand_~a {~%" num)
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
    (format dest "}~%")))


#| testing code

	(defparameter result (parse-json *j*))

	result

	(/= -1 -1 -2 -3 0)


	(scaffold-vec (gethash 0 *origami*))
	|#

