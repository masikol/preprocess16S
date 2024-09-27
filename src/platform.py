# -*- coding: utf-8 -*-
# Module contains stuff for behaving differently on different platforms.

import sys

def platf_depend_exit(exit_code):
    # Function asks to press ENTER press on Windows
    #     and exits after that.
    # preprocess16S cannot now be executed on Windows, but anyway.
    # :type exit_code: int;
    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit
