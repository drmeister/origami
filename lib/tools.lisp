(in-package :cl-jupyter-user)

;;; Start a swank server
(defun slime ()
  (load "/home/app/slime/swank-loader.lisp")
  (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
    (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
  (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
    (mp:process-run-function 'swank-main
			     (lambda () (funcall swank-create-server
						 :port 4005
						 :interface "0.0.0.0")))))

