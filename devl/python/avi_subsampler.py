from pathlib import Path
import math

import cv2
from tqdm import tqdm

INPUT_PATH = Path("C:/Users/dirk_claudia/Pictures/dji/igel_2026/DJI_20261030152457_0005_D.MP4")
OUTPUT_PATH = Path("C:/Users/dirk_claudia/Pictures/dji/igel_2026/DJI_20261030152457_0005_D.20.avi")
FRAME_INTERVAL = 20


def main() -> None:
    if not INPUT_PATH.is_file():
        raise FileNotFoundError(f"Video not found: {INPUT_PATH}")
    if INPUT_PATH.resolve() == OUTPUT_PATH.resolve():
        raise ValueError("Input and output paths must be different.")
    if FRAME_INTERVAL < 1:
        raise ValueError("FRAME_INTERVAL must be at least 1.")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    capture = cv2.VideoCapture(str(INPUT_PATH))
    if not capture.isOpened():
        raise RuntimeError(f"Could not open video: {INPUT_PATH}")

    writer = None
    try:
        input_fps = capture.get(cv2.CAP_PROP_FPS)
        reported_frame_count = int(capture.get(cv2.CAP_PROP_FRAME_COUNT))
        width = int(capture.get(cv2.CAP_PROP_FRAME_WIDTH))
        height = int(capture.get(cv2.CAP_PROP_FRAME_HEIGHT))
        if not math.isfinite(input_fps) or input_fps <= 0 or width <= 0 or height <= 0:
            raise RuntimeError("The input video has invalid frame-rate or dimensions.")

        output_fps = input_fps
        writer = cv2.VideoWriter(
            str(OUTPUT_PATH),
            cv2.VideoWriter_fourcc(*"MJPG"),
            output_fps,
            (width, height),
        )
        if not writer.isOpened():
            raise RuntimeError(f"Could not create AVI file: {OUTPUT_PATH}")

        frames_read = 0
        frames_written = 0
        with tqdm(
            total=int(reported_frame_count / FRAME_INTERVAL) or None,
            desc="Verarbeite Frames",
            unit="Frame",
        ) as progress:
            while True:
                # capture.set(cv2.CAP_PROP_POS_FRAMES, frames_read)

                for _ in range(FRAME_INTERVAL):
                    success = capture.grab()
                    if not success:
                        break
                success, frame = capture.retrieve()
                if not success:
                    raise RuntimeError(f"Could not retrieve frame {frames_read}.")
                writer.write(frame)
                frames_written += 1
                frames_read += FRAME_INTERVAL
                progress.update(1)

        print(f"Expected {reported_frame_count} frames, but read {frames_read}.")
    finally:
        print("finally error occurs.")

    capture.release()
    if writer is not None:
        writer.release()

    print(f"Frames read: {frames_read}")
    print(f"Frames written: {frames_written}")
    print(f"Output frame rate: {output_fps:.6f} fps")
    print(f"Approximate output duration: {frames_written / output_fps:.2f} seconds")
    print(f"Output: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
