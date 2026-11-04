"""Development-only, explicit scientific decision; no score optimizer."""
import experiment as e

REASON = (
    'Retain five-term B32. A failed its targeted low-M repair and worsened primary scores; B lost substantial predictive information. '
    'C has the expected negative selectivity coefficient in all 47 LOCO fits and improves loss in 32/47 chemistries (median delta -0.00659), '
    'which makes the hypothesis scientifically interesting. However, its pooled AP gain is only 0.001386 (paired chemistry-bootstrap 95% interval [-0.02208, 0.01812]); '
    'LOCO loss worsens by 0.001618 (interval [-0.01395, 0.01839]), macro chemistry loss worsens from 0.38199 to 0.39448, '
    'and complementary grouped-CV loss worsens in 16/20 repeats. The losses for AcDOT-[OE3]2, D18, and several DOT chemistries outweigh many small gains. '
    'Training improvement and a stable coefficient sign do not establish transfer benefit. These small and inconsistent gains do not justify a sixth term. '
    'This is a negative equation-repair result, not evidence that B32 is a physical law. All three allowed new structures have been used; '
    'no new candidates or structure decisions will follow the holdout.'
)

if __name__ == '__main__':
    e.log('Final development feedback and selection reasoning', REASON)
    e.freeze('B32', REASON)
