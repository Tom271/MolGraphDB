from src.hello import main


def test_main(capsys):
    main()
    captured = capsys.readouterr()
    assert captured.out == "Hello from python-project-template!\n"
