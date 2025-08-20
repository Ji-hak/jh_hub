
from model import MarketModel

def main():
    """
    시뮬레이션을 설정하고 실행하는 메인 함수입니다.
    """
    # 시뮬레이션 파라미터 설정
    num_households = 50  # 가계 에이전트 수
    num_firms = 10       # 기업 에이전트 수
    simulation_rounds = 50 # 총 라운드 수

    # MarketModel 객체 생성
    simulation = MarketModel(
        num_households=num_households,
        num_firms=num_firms,
        rounds=simulation_rounds
    )

    # 시뮬레이션 실행
    simulation.run()

    # 시뮬레이션 종료 및 결과 처리
    simulation.finish()


if __name__ == '__main__':
    # 이 스크립트가 직접 실행될 때 main() 함수를 호출합니다.
    main()
