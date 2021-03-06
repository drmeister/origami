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
(defun parse-json (json csv-data)
  (format t "Starting parse-json~%")
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
    ;;deal with sequence
    (format t "here~%")
    (unless csv-data
      (warn "No CSV; nucleobase information will be defaulted"))
    (format t "About to assign node-staple-sequences~%")
    (let ((node-staple-sequences (make-hash-table :test #'eq)))
      (loop for (num pos sequence) in csv-data
	 for start-node = (lookup-node vstrands num pos :staple-vec)
	 do (setf (gethash start-node node-staple-sequences) sequence))
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
      ;;fill base names
      (format t "About to apply node-staple-sequences~%")
      (loop for sequence being the hash-values in node-staple-sequences using (hash-key start-node)
	 do (loop with node = start-node
	       for base across sequence
	       do (unless (char= base #\?)
		    (setf (base-name node) base))
		 (setf node (forward-node node))
		 (unless node (return))))
      vstrands)))

   
(defun read-csv-sequence-file (csv-filename)
  (if (probe-file csv-filename)
      (with-open-file (csv-filename :direction :input)
	(read-csv:parse-csv csv-file))
      nil))

(defun find-title-row (data)
  (find-if (lambda (row)
	     (alpha-char-p (char (first row) 0)))
	   data))

(defun title-positions (title-row)
  (values (position "Start" title-row :test #'string=)
	  (position "Sequence" title-row :test #'string=)))

(defun positions-and-sequences (csv)
  (let* ((title-row (find-title-row csv))
	 (data (remove title-row csv)))
    (multiple-value-bind (start sequence)
	(title-positions title-row)
      (loop for row in data
	 collect (multiple-value-bind (strand pos)
		     (position-from-string (nth start row))
		   (list strand pos
			 (nth sequence row)))))))

(defun position-from-string (string)
  (multiple-value-bind (first stop)
      (parse-integer string :junk-allowed t)
    (values first (parse-integer string :start (1+ stop) :junk-allowed t))))

(defun apply-default-sequence-to-cstrands (cstrands)
  (loop for strand in cstrands
     do (loop for node in (dna-strand-as-list-of-nodes strand)
	   do (progn
		(setf (base-name node) :G)
		(when (hbond-node node)
		  (setf (base-name (hbond-node node)) :C))))))

(defun node-has-identifier-p (node)
  (slot-boundp node 'identifier))

(defun node-has-base-name-p (node)
  (slot-boundp node 'base-name))

(defun complementary-base (base)
  (ecase base
    ((:c) :g)
    ((:g) :c)
    ((:a) :t)
    ((:t) :a)))

(defun fill-empty-bases-with-default (cstrands)
  ;;; Apply the sequence information in csv-sequence to the vstrands
  ;;;    If there is no csv-sequence then use :G/:C
  (loop for strand in cstrands
     do (loop for node in (dna-strand-as-list-of-nodes strand)
	   for hnode = (hbond-node node)
	   do (format t "fill-empty-bases-with-default node -> ~a~%" node)
	   do (if (node-has-base-name-p node)
		  (progn
		    (format t "TRUE node-has-base-name-p~%")
		    (cond
		      ((not hnode))
		      ((node-has-base-name-p hnode)
		       (unless (eql (base-name hnode)
				    (complementary-base (base-name node)))
			 (error "Base pair doesn't match: ~a" node)))
		      (t (setf (base-name hnode)
			       (complementary-base (base-name node))))))
		  (progn
		    (format t "NOT node-has-base-name-p~%")
		    (cond ((not hnode)
			 (setf (base-name node) :G))
			((node-has-base-name-p hnode)
			 (setf (base-name node)
			       (complementary-base (base-name hnode))))
			(t (setf (base-name node) :G
				 (base-name hnode) :C))))))))

(defun parse-cadnano (file-name)
  (let ((json-filename (make-pathname :type "json" :defaults file-name))
	(csv-filename (make-pathname :type "csv" :defaults file-name)))
    (let* ((json (with-open-file (json-file json-filename :direction :input)
		   (json:decode-json json-file)))
	   (csv (with-open-file (csv-file csv-filename :direction :input
					  :if-does-not-exist nil)
		  (positions-and-sequences
		   (read-csv:parse-csv csv-file))))
	   (vstrands (parse-json json csv))
	   (cstrands (single-or-double-classified-strands vstrands)))
      (format t "About to fill-empty-bases-with-default~%")
      (fill-empty-bases-with-default cstrands)
      cstrands)))
    
(defun BUILD-NODE (strand-json vec num strand-name)
  (loop for index from 0 below (length vec)
     for node-json in strand-json
     unless (apply #'= -1 node-json)
     do (setf (elt vec index) (make-instance 'node :identifier (list :j num strand-name index)))))

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
   (base-name :initarg :base-name :accessor base-name)
   (residue :initform nil :initarg :residue :accessor residue)
   (strand :initarg :strand :accessor strand)
   (identifier :initarg :identifier :reader identifier)))

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
	    (setf (elt new-vec dest-cursor) (make-instance 'node :identifier (list :A num strand-name s-l-cursor :I (elt loop-vec s-l-cursor))))
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
    (format dest "\"~a\" -> \"~a\" [color=~a];~%" (identifier node) (identifier next-node) color)))

(defun draw-connection (dest node)
  (safe-draw-link dest node (forward-node node) "violet")
  (safe-draw-link dest node (backward-node node) "blue")
  (safe-draw-link dest node (hbond-node node) "orange"))

(defun draw-node (dest node)
;  (let ((*print-pretty* nil)))
  (format dest "\"(~{~a~^ ~})\"->" (identifier node)))

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

(defun load-*bases* (dir)
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
	    atom	      ; treat everything like a carbon for now
	    1.908	      ; CT from parm99.dat
	    0.1094	      ; CT from parm99.dat
	    0.0		      ; no charge
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
	      atom	      ; treat everything like a carbon for now
	      1.908	      ; CT from parm99.dat
	      0.1094	      ; CT from parm99.dat
	      0.0	      ; no charge
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
  
