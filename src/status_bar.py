# -*- coding: utf-8 -*-

def status_bar(iteration_func, infpaths):

    nreads = src.fastq.count_reads(infpaths)
    bar_len = int(os.get_terminal_size().columns * 0.40)
    next_print_num = int(nreads * 0.01)
    inc_num = next_print_num
    i = 0
    sys.stdout.write(' [>{}] 0/{} (0%)'.format(' '*bar_len, nreads))

    for fastq_records in src.fastq.fastq_generator(infpaths):
        iteration_func(fastq_records)

        i += 1

        # Update status bar
        if i > next_print_num:
            bar_len = int(os.get_terminal_size().columns * 0.40)
            done_ratio = i / nreads
            sys.stdout.write('\r [{}>{}] {}/{} ({}%)'.format('='*int(bar_len*done_ratio),
                ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)))
            next_print_num += inc_num
        # end if
    # end for

    # Update status bar
    bar_len = int(os.get_terminal_size().columns * 0.40)
    done_ratio = i / nreads
    sys.stdout.write('\r [{}>{}] {}/{} ({}%)\n'.format('='*int(bar_len*done_ratio),
        ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)))
    next_print_num += inc_num
# end def status_bar
