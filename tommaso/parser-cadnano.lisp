
(ql:quickload "cl-json")


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
       do (BUILD-NODE staple-json (staple-vec vstrand))	 
					;       do (BUILD-NODE (scaffold-vec vstrand) scaffold-json)
       do (BUILD-NODE scaffold-json (scaffold-vec vstrand))
       do (intra-helix-connect (staple-vec vstrand) (scaffold-vec vstrand)))
    ;; connect the nodes
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for staple-json = (staple-json vstrand)
       for staple-vec = (staple-vec vstrand)
       for scaffold-json = (scaffold-json vstrand)
       for scaffold-vec = (scaffold-vec vstrand)
       do (CONNECT-EVERYTHING vstrands num staple-json staple-vec #'staple-vec)
       do (CONNECT-EVERYTHING vstrands num scaffold-json scaffold-vec #'scaffold-vec))
    ;; Deal with loop
    vstrands))
    
(defun BUILD-NODE (strand-json vec)
  (loop for index from 0 below (length vec)
     for node-json in strand-json
     unless (apply #'= -1 node-json)
     do (setf (elt vec index) (make-instance 'node :name (cons :json index)))))

(defun INTRA-HELIX-CONNECT (staple-vec scaffold-vec)
  (loop for index from 0 below (length staple-vec)
     for staple-node = (elt staple-vec index)
     for scaffold-node = (elt scaffold-vec index)
     do (when (and staple-node scaffold-node)
	  (setf (hbond-node staple-node) scaffold-node)
	  (setf (hbond-node scaffold-node) staple-node))))
    
(defun CONNECT-EVERYTHING (vstrands num strand-json strand-vec accessor)
  (loop for json-info in strand-json
     for index from 0 below (length strand-json)
     for node = (elt strand-vec index)
     when node
     do (let* ((fwd-vstrand (first json-info))
	       (fwd-pos     (second json-info))
	       (bkd-vstrand (third json-info))
	       (bkd-pos     (fourth json-info)))
	  (when (/= fwd-vstrand -1)
	    (let ((forward-node (LOOKUP-NODE vstrands fwd-vstrand fwd-pos accessor)))
	      (assert forward-node)
	      (setf (forward-node node) forward-node)))
	  (when (/= bkd-vstrand -1)
	    (let ((back-node (LOOKUP-NODE vstrands bkd-vstrand bkd-pos accessor)))
	      (assert back-node)
	      (setf (back-node node) back-node))))))
	      
(defclass node ()
  ((hbond-node :initarg :hbond-node :accessor hbond-node)
   (forward-node :initarg :forward-node :accessor forward-node)
   (back-node :initarg :back-node :accessor back-node)
   (name :initarg :name :reader name)))

(defun LOOKUP-NODE (vstrands vstrand-num pos accessor)
  (let* ((vstrand (gethash vstrand-num vstrands))
	 (strand (funcall accessor vstrand)))
    (elt strand pos)))

(defun skip-procedure (x)
  (let ((xf (forward x))
	(xb (backward x)))
    ((setf (forward xb) xf)
     (setf (backward xf) xb)
     (setf x NIL))))
(defun arrow-direction (vec)
  (let ((step-direction (loop for index from 0 below (length vec)
			   for node = (elt vec index)
			   when node
			   collect (let (fwd (forward node))
				     (when fwd
				       (- (position fwd vec) index))))))
    (cond
      (((every (lambda (x) (if x (plusp x) T)) step-direction))
       1)
      (((every (lambda (x) (if x (minusp x) T)) step-direction))
       -1)
      (t (error "arrow direction in the vector is neither foward nor backward")))
    ))
				  
					     
(defun skip-loop (skip loop staple-vec scaffold-vec)
  ;;must check if 4 input have the same lenght
  (let* ((old-vec-min-length
	  ;;in presence of gaps...
					;	  (MAX (loop for node across staple-vec
					;				when node
					;		  count)
					;	       (loop for node across scaffold-vec
					;				when node
					;		  count))
	  (loop for index from 0 below (length scaffol-vec)
	     when (or (elt scaffold-vec index)  (elt staple-vec index))
	     count))
	 (new-vec-length (apply #'+ old-vec-min-length (apply #'+ skip) loop))
	 (new-scaf-vec (make-array new-vec-length))
	 (new-stap-vec (make-array new-vec-length))
	 (s-l-cursor 0)
	 (dest-cursor 0)
	 (stap-direction (arrow-direction staple-vec))
	 (scaf-direction (arrow-direction scaffold-vec)))
	       (cond ((and (= stap-direction 1)
			 (= scaf-direction -1))
					;stap->f
					;scaf->b
		       )
		     ((and (= stap-direction -1)
			   (= scaf-direction 1))
					;stap->b
					;scaf->f
		      )

    (loop with s-l-cursor = 0
       with dest-cursor = 0
	 with length-skip = (length skip)
       until (= s-l-cursor length-skip)
	 do (cond
	      ((and (= (elt skip s-l-cursor) -1)
		    (= (elt loop s-l-cursor) 0))
	       ;;do a skip procedure
	       (skip-procedure (elt staple-vec s-l-cursor))
	       (skip-procedure (elt scaffold-vec s-l-cursor))
	       (incf s-l-cursor)
	       )
	      ((and (> (elt skip s-l-cursor) 0)
		    (= (elt loop s-l-cursor) 0))
	       ;;do an insertion
			 ;;forward arrow		    
	       (let ((New-node (elt new-f-vec dest-cursor))
		     (Doubled-node (elt f-vec s-l-cursor))
		     (b-d-node (backward (elt f-vec s-l-cursor))))
		 (setf New-node (make-instance 'node :name (cons :loop (elt loop s-l-cursor))))
		 (setf (backward Doubled-node) New-node)
		 (setf (forward New-node) Doubled-node)
		 (setf (backward New-node) b-d-node)
		 (setf (forward b-d-node) New-node)
		 )
	       ;;backward arrow
	       (let ((New-node (elt new-b-vec dest-cursor))
		     (Doubled-node (elt b-vec s-l-cursor))
		     (f-d-node (forward (elt b-vec s-l-cursor))))
		 (setf New-node (make-instance 'node :name (cons :loop (elt loop s-l-cursor))))
		 (setf (forward Doubled-node) New-node)
		 (setf (backward New-node) Doubled-node)		 
		 (setf (forward New-node) f-d-node)
		 (setf (backward b-d-node) New-node)
		 )
	       (decf (elt loop s-l-cursor))
	       (incf dest-cursor)
	       )
	      ((and (= (elt skip s-l-cursor) 0)
		(= (elt loop s-l-cursor) 0))
	       ;;do a copy
	       (incf dest-cursor)
	       (incf s-l-cursor)
	       )
	      (t (error "An impossible step was encountered")))
  
#| testing code

(defparameter result (parse-json *j*))

result

(/= -1 -1 -2 -3 0)

|#
